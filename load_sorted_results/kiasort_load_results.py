"""Load kiaSort outputs into Python."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import h5py
import scipy.io as sio


def kiasort_load_results(
    results_path,
    curated: bool = False,
    individual_waveform: bool = False,
) -> Dict[str, Any]:
    """Load kiaSort outputs.

    Parameters
    ----------
    results_path : str | Path
        The outputFolder passed to the sorter.
    curated : bool
        Prefer RES_Sorted/<field>_curated.h5; fall back to non-curated
        with a warning if missing.
    individual_waveform : bool
        Also load per-spike waveforms when saved.

    Returns
    -------
    dict with keys: source, spike_idx, channelNum, unifiedLabels,
    [labelOnChannel], [waveforms], units (list of dicts with id,
    channel, meanWaveform, isolation).
    """
    results_path = Path(results_path)
    res_folder = results_path / "RES_Sorted"
    if not res_folder.is_dir():
        raise FileNotFoundError(f"RES_Sorted not found at {res_folder}")

    suffix = ""
    if curated:
        if (res_folder / "spike_idx_curated.h5").exists() and \
           (res_folder / "unifiedLabels_curated.h5").exists():
            suffix = "_curated"
        else:
            warnings.warn(
                f"Curated outputs missing in {res_folder}. Loading non-curated."
            )

    out: Dict[str, Any] = {
        "source":        "RES_Sorted" + suffix,
        "spike_idx":     _read_h5(res_folder, f"spike_idx{suffix}"),
        "channelNum":    _read_h5(res_folder, f"channelNum{suffix}"),
        "unifiedLabels": _read_h5(res_folder, f"unifiedLabels{suffix}"),
    }

    labels_file = res_folder / "labels.h5"
    if labels_file.exists():
        out["labelOnChannel"] = _read_h5(res_folder, "labels")

    if individual_waveform:
        wf_file = res_folder / f"waveforms{suffix}.h5"
        if wf_file.exists():
            with h5py.File(wf_file, "r") as f:
                out["waveforms"] = np.asarray(f[f"/waveforms{suffix}"])
        else:
            warnings.warn(f"Waveform file {wf_file} not found.")

    out["units"] = _load_units(results_path, res_folder, suffix)
    return out


def _read_h5(folder: Path, name: str) -> np.ndarray:
    with h5py.File(folder / f"{name}.h5", "r") as f:
        return np.asarray(f[f"/{name}"]).flatten()


def _load_units(results_path: Path, res_folder: Path, suffix: str) -> List[Dict[str, Any]]:
    units: List[Dict[str, Any]] = []

    curated_mat = res_folder / "curated_sample.mat"
    sorted_samp = results_path / "Sorted_Samples" / "sorted_samples.mat"

    if suffix == "_curated" and curated_mat.exists():
        d = sio.loadmat(str(curated_mat), simplify_cells=True)
        c = d.get("curatedSamples")
        if isinstance(c, dict) and "unifiedLabels" in c:
            ids   = np.atleast_1d(c["unifiedLabels"]).flatten()
            chans = np.atleast_1d(c.get("channelNum", np.full(ids.shape, np.nan))).flatten()
            waves = c.get("waveform")
            iso   = c.get("unitIsolation")
            for k, lab in enumerate(ids):
                u: Dict[str, Any] = {"id": _to_py(lab)}
                if k < chans.size:
                    u["channel"] = _to_py(chans[k])
                if isinstance(waves, np.ndarray) and waves.ndim >= 3 and k < waves.shape[0]:
                    u["meanWaveform"] = np.asarray(waves[k])
                if iso is not None:
                    iso_arr = np.atleast_1d(iso).flatten()
                    if k < iso_arr.size:
                        u["isolation"] = str(iso_arr[k])
                units.append(u)
            return units

    if sorted_samp.exists():
        d = sio.loadmat(str(sorted_samp), simplify_cells=True)
        ccs = d.get("crossChannelStats")
        if isinstance(ccs, dict) and "unified_labels" in ccs:
            uL = ccs["unified_labels"]
            ids   = np.atleast_1d(uL["label"]).flatten()
            chans = np.atleast_1d(uL.get("channelID", np.full(ids.shape, np.nan))).flatten()
            waves = uL.get("meanWaveforms")
            for k, lab in enumerate(ids):
                u = {"id": _to_py(lab)}
                if k < chans.size:
                    u["channel"] = _to_py(chans[k])
                if isinstance(waves, np.ndarray) and waves.ndim >= 3 and k < waves.shape[0]:
                    u["meanWaveform"] = np.asarray(waves[k])
                units.append(u)

    return units


def _to_py(x):
    try:
        if np.isnan(x):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(x, "item"):
        return x.item()
    return x
