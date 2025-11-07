"""
PyRMS example session: calculating RMS statistics for a multi-state object in PyMOL 
1A1T: https://doi.org/10.1126/science.279.5349.384
"""

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
