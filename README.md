# PyRMS

PyRMS is a small Python utility and PyMOL plugin for calculating pairwise RMS deviations across states of multistate structures with detailed reporting for both all atoms and heavy atoms (hydrogens excluded). It supports interactive use inside PyMOL and scripted usage from Python.

## Features

- Compute RMS deviations per state and for arbitrary atom selections\
- Automatic calculation for both "all" atoms and "heavy" atoms (no hydrogens)\
- Per-state and per-pair RMS reporting (summary and full detail modes)\
- Optional export of timestamped text reports\
- Configurable verbosity and display levels

## Online use (no install)

The script can be modified and executed [there](https://mybinder.org/v2/gh/EricLarG4/pyRMS/977364510360a990ec76e95537c77a35a109eaad?urlpath=lab%2Ftree%2Fpairwise_rms.ipynb).

## Installation

Choose one of the following install workflows.

### For PyMOL GUI use

1. Open PyMOL (version 2.0 or higher recommended).
2. Download `rms_calculator.py` and place it in a directory accessible by PyMOL.
3. Load the script from within PyMOL (File → Open).
4. The function `calculate_rms_stats` can then be called from the PyMOL command line.

Alternatively, you can proceed without downloading the script, and replace steps
2 and 3 by pasting the following command directly into the PyMOL command line:

```python
run https://raw.githubusercontent.com/EricLarG4/pyRMS/refs/heads/main/rms_calculator.py
```

### For scripted/IDE use

1. Ensure you have Python 3.11+ installed.
2. Download or clone the repository.
3. Create an environment and install dependencies *via*:

Conda (recommended for standalone/scripted use)

```bash
conda env create -f environment.yml
conda activate pyrms
```

There is also a pyproject.toml for pip/poetry users, however PyMOL installation requires conda.

## Parameters

| Parameter         | Type | Default   | Description                                                        |
| ----------------- | ---- | --------- | ------------------------------------------------------------------ |
| `obj_to_fit`    | str  | —        | PyMOL object name or selection to analyze (required)               |
| `selection`     | str  | `"all"` | Additional PyMOL selection string                                  |
| `quiet`         | bool | `False` | Suppress printing when `True`                                    |
| `display_level` | int  | `1`     | 0 → summary only, 1 → per-state means, 2 → full per-pair detail |
| `export_report` | bool | `False` | Export timestamped report to a file when `True`                  |
| `report_path`   | str  | `.`     | Directory to save exported report                                  |

## Returns

The function returns a dictionary with RMS statistics and metadata. Typical keys:

- `modes`: detailed stats for `"all"` and `"no_hydrogen"` modes\
- `selection`: resolved selection string used for calculations\
- `atom_counts`: atom counts per mode\
- `states`: number of states analyzed\
- `report_file`: path to the exported report (if `export_report=True`)\
- `error`: optional error messages (present only on failure)

Example access:

```python
mean_all = result["modes"]["all"]["overall_mean"]
state_means = result["modes"]["no_hydrogen"].get("state_means")
```

## Usage

Call the `calculate_rms_stats` function from within PyMOL or a Python script.

If called without arguments, `calculate_rms_stats()` will display help information.

### Interactive (PyMOL GUI)

The minimal usage from the PyMOL command line is as follows:

```python
calculate_rms_stats object_name
```

Where `object_name` is the name of the multistate object loaded in PyMOL.
This will output RMS statistics for all atoms and heavy atoms (no hydrogens)
with default settings.
**No chain/residue/atom filtering** is applied (selection = "all"), and no report file is generated.

Below is an example of advanced use for the 1A1T structure:

```python
# Load example structure
cmd.fetch("1A1T")

# Calculate RMS stats for all atoms (+/- hydrogens)
full_structure = calculate_rms_stats(\
    "1A1T",\
    selection="all",\
    quiet=False,\
    display_level=1,\
    export_report=False,\
    report_path="."\
)

# Calculate RMS stats for RNA only (chain A)
rna_only = calculate_rms_stats(\
    "1A1T",\
    selection="resn A+G+C+U",\
    quiet=False,\
    display_level=1,\
    export_report=False,\
    report_path="."\
)

# Caculate RMS stats for protein only (chain B)
protein_only = calculate_rms_stats(\
    "1A1T",\
    selection="chain A and not resn ZN",\
    quiet=False,\
    display_level=1,\
    export_report=False,\
    report_path="."\
)

# Calculate RMS stats for protein backbone only (chain B + backbone atoms)
protein_backbone = calculate_rms_stats(\
    "1A1T",\
    selection="chain A and backbone",\
    quiet=False,\
    display_level=1,\
    export_report=False,\
    report_path="."\
)

# Access summary values

## Full structure overall mean RMS
mean_all = full_structure["modes"]["all"]["overall_mean"]
sd_all = full_structure["modes"]["all"]["overall_sd"]
print(f"Overall mean RMS (all atoms): {mean_all:.2f} +/- {sd_all:.2f} Å")

## Full structure overall mean RMS (no hydrogens)
mean_no_hydrogen = full_structure["modes"]["no_hydrogen"]["overall_mean"]
sd_no_hydrogen = full_structure["modes"]["no_hydrogen"]["overall_sd"]
print(f"Overall mean RMS (no hydrogens): {mean_no_hydrogen:.2f} +/- {sd_no_hydrogen:.2f} Å")

# Compare mean RMS between state 1 and state 10 for protein only
state_1_to_10_rms = protein_only["modes"]["all"]["per_pair_rms"][(1, 10)]
print(f"RMS between state 1 and state 10 (protein only): {state_1_to_10_rms:.2f} Å")
```

### IDE

Examples from a Python script, available in the `examples/` directory. Reports will be saved to `./examples/reports/`.

```python
# Import necessary modules
from pymol import cmd
from rms_calculator import calculate_rms_stats as calc_rms

# Load example structure
cmd.reinitialize("everything") # Clear existing data to avoid conflicts (PyMOL error: loading mmCIF into existing object not supported)
cmd.fetch("1A1T")

# Calculate RMS stats for all atoms (+/- hydrogens)
full_structure = calc_rms(
    "1A1T",
    selection="all",
    quiet=False,
    display_level=2,
    export_report=True,
    report_path="./examples/reports"
)

# Calculate RMS stats for RNA only (chain A)
rna_only = calc_rms(
    "1A1T",
    selection="resn A+G+C+U",
    quiet=False,
    display_level=2,
    export_report=True,
    report_path="./examples/reports"
)

# Caculate RMS stats for protein only (chain B)
protein_only = calc_rms(
    "1A1T",
    selection="chain A and not resn ZN",
    quiet=False,
    display_level=2,
    export_report=True,
    report_path="./examples/reports"
)

protein_backbone = calc_rms(
    "1A1T",
    selection="chain A and backbone",
    quiet=False,
    display_level=2,
    export_report=True,
    report_path="./examples/reports"
)

# Access summary values

## Full structure overall mean RMS
mean_all = full_structure["modes"]["all"]["overall_mean"]
sd_all = full_structure["modes"]["all"]["overall_sd"]
print(f"Overall mean RMS (all atoms): {mean_all:.2f} +/- {sd_all:.2f} Å")

## Full structure overall mean RMS (no hydrogens)
mean_no_hydrogen = full_structure["modes"]["no_hydrogen"]["overall_mean"]
sd_no_hydrogen = full_structure["modes"]["no_hydrogen"]["overall_sd"]
print(f"Overall mean RMS (no hydrogens): {mean_no_hydrogen:.2f} +/- {sd_no_hydrogen:.2f} Å")

```

## Report output

- When `export_report=True`, PyRMS writes a timestamped text report to `report_path` (default: current working directory).\
- The report contains per-state summaries and, if `display_level=2`, full per-pair RMS details. The returned dictionary includes `report_file` with the report path when exported.

## License

MIT License — feel free to use, modify, and redistribute under the terms of the MIT license.
