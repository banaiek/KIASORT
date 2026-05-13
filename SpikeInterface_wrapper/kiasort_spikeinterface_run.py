"""Example entry point for running kiaSort on a SpikeInterface recording.

Edit the paths, recording loader, and config_overrides below, then run.
"""

from pathlib import Path

import spikeinterface.extractors as se

from kiasort_spikeinterface import run_kiasort


# ---- Paths --------------------------------------------------------------
kiasort_path  = "/path/to/Nogui_Files"
output_folder = "/path/to/results"
matlab_bin    = "matlab"

# ---- Recording ----------------------------------------------------------
# Replace with whichever SI extractor matches your data. Examples:
#   recording = se.read_binary("rec.dat", num_channels=128,
#                              sampling_frequency=30000, dtype="int16")
#   recording = se.read_spikeglx("/path/to/spikeglx_folder")
#   recording = se.read_openephys("/path/to/openephys_folder")
recording = se.read_binary(
    file_paths="recording.dat",
    sampling_frequency=30000,
    num_channels=128,
    dtype="int16",
)
# If the extractor does not set channel locations, attach a probe:
# from probeinterface import generate_linear_probe
# recording = recording.set_probe(generate_linear_probe(num_elec=128, ypitch=20))

# ---- Config overrides (same fields as kiaSort_main_noGUI_run.m) --------
config_overrides = {
    "bandpass":             [500, 6000],
    "commonRef":            "median",         # 'median' | 'mean' | 'none'
    "denoising":            True,
    "extremeNoise":         False,
    "min_threshold":        5.5,
    "waveform_radius":      150,
    "num_channel_extract":  5,
    "numSampleChunks":      300,
    "sortingChunkDuration": 120,
    "sort_only":            False,
    "useGPU":               True,
    "parallelProcessing":   False,
    "extractWaveform":      False,
}

# ---- Run ----------------------------------------------------------------
sorting = run_kiasort(
    recording=recording,
    output_folder=output_folder,
    kiasort_path=kiasort_path,
    matlab_bin=matlab_bin,
    config_overrides=config_overrides,
    keep_intermediate=False,
    verbose=True,
)

print(f"Done. {len(sorting.unit_ids)} units, "
      f"{sum(len(sorting.get_unit_spike_train(u)) for u in sorting.unit_ids)} spikes.")
print(f"Results: {Path(output_folder) / 'RES_Sorted'}")
