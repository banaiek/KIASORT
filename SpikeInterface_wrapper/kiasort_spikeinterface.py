"""Run kiaSort on a SpikeInterface recording.

Requires MATLAB on PATH (or via `matlab_bin`) and the kiaSort Nogui_Files
directory available on disk.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional, Mapping, Any

import numpy as np
import h5py
from scipy.io import savemat

import spikeinterface as si
from spikeinterface.core import write_binary_recording


def run_kiasort(
    recording,
    output_folder,
    kiasort_path,
    matlab_bin: str = "matlab",
    config_overrides: Optional[Mapping[str, Any]] = None,
    keep_intermediate: bool = False,
    verbose: bool = True,
):
    """Run kiaSort (no-GUI) on a SpikeInterface recording.

    Parameters
    ----------
    recording : si.BaseRecording
        Probe-equipped recording.
    output_folder : str | Path
        Where the .dat, channel map, and kiaSort outputs go.
    kiasort_path : str | Path
        Path to the kiaSort Nogui_Files directory.
    matlab_bin : str
        MATLAB launcher command.
    config_overrides : dict, optional
        Fields applied on top of the defaults inside run_kiasort_nogui.
        samplingFrequency / numChannels / dataType are filled in from
        the recording if not provided here.
    keep_intermediate : bool
        Keep the intermediate .dat and channel map after the run.
    verbose : bool
        Echo MATLAB output.

    Returns
    -------
    sorting : si.NumpySorting
    """
    output_folder = Path(output_folder).resolve()
    output_folder.mkdir(parents=True, exist_ok=True)
    kiasort_path = Path(kiasort_path).resolve()

    fs = float(recording.get_sampling_frequency())
    n_ch = int(recording.get_num_channels())
    try:
        locs = np.asarray(recording.get_channel_locations(), dtype=float)
    except Exception:
        locs = np.zeros((n_ch, 2), dtype=float)
        locs[:, 1] = np.arange(n_ch)
    if locs.ndim == 1:
        locs = locs[:, None]
    if locs.shape[1] < 2:
        locs = np.column_stack([locs, np.arange(n_ch)])

    dat_path = output_folder / "recording.dat"
    write_binary_recording(
        recording, file_paths=[str(dat_path)], dtype="int16", verbose=False
    )

    chan_map_path = output_folder / "channel_map.mat"
    savemat(str(chan_map_path), {
        "chanMap":    np.arange(1, n_ch + 1, dtype=np.int32),
        "connected":  np.ones(n_ch, dtype=bool),
        "xcoords":    locs[:, 0].astype(float),
        "ycoords":    locs[:, 1].astype(float),
    })

    overrides = dict(config_overrides or {})
    overrides.setdefault("samplingFrequency", fs)
    overrides.setdefault("numChannels", n_ch)
    overrides.setdefault("dataType", "int16")

    script_path = output_folder / "_run_kiasort.m"
    script_path.write_text(_build_matlab_script(
        kiasort_path, dat_path, output_folder, chan_map_path, overrides
    ))

    cmd = [matlab_bin, "-batch", f"run('{script_path.as_posix()}')"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if verbose:
        if res.stdout:
            print(res.stdout)
        if res.stderr:
            print(res.stderr)
    if res.returncode != 0:
        raise RuntimeError(
            f"MATLAB exited with code {res.returncode}\n{res.stderr}"
        )

    res_dir = output_folder / "RES_Sorted"
    spike_idx_path = res_dir / "spike_idx.h5"
    labels_path = res_dir / "unifiedLabels.h5"
    if not spike_idx_path.exists() or not labels_path.exists():
        raise FileNotFoundError(
            f"kiaSort outputs not found in {res_dir}"
        )

    with h5py.File(spike_idx_path, "r") as f:
        spike_idx = np.asarray(f["/spike_idx"]).flatten().astype(np.int64)
    with h5py.File(labels_path, "r") as f:
        labels = np.asarray(f["/unifiedLabels"]).flatten().astype(np.int64)

    keep = labels >= 0
    spike_idx = spike_idx[keep]
    labels = labels[keep]

    sorting = si.NumpySorting.from_times_labels(
        times_list=spike_idx, labels_list=labels, sampling_frequency=fs
    )

    if not keep_intermediate:
        for p in (dat_path, chan_map_path, script_path):
            try:
                p.unlink()
            except FileNotFoundError:
                pass

    return sorting


def _format_matlab_value(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        s = v.replace("'", "''")
        return f"'{s}'"
    if isinstance(v, (list, tuple, np.ndarray)):
        arr = np.asarray(v).ravel()
        parts = [_format_matlab_value(x.item() if hasattr(x, "item") else x) for x in arr]
        return "[" + " ".join(parts) + "]"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    return repr(v)


def _build_matlab_script(kiasort_path: Path, dat_path: Path,
                        output_folder: Path, chan_map_path: Path,
                        overrides: Mapping[str, Any]) -> str:
    lines = [
        f"addpath(genpath('{Path(kiasort_path).as_posix()}'));",
        "cfg_overrides = struct();",
    ]
    for k, v in overrides.items():
        lines.append(f"cfg_overrides.{k} = {_format_matlab_value(v)};")
    lines.append(
        f"run_kiasort_nogui('{dat_path.as_posix()}', "
        f"'{output_folder.as_posix()}', "
        f"'{chan_map_path.as_posix()}', cfg_overrides);"
    )
    return "\n".join(lines) + "\n"
