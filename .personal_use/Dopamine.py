"""
PyRMS dopamine session
"""

# Import necessary modules
from pymol import cmd
from rms_calculator import calculate_rms_stats as calc_rms

# Load example structure
cmd.reinitialize # Clear existing data to avoid conflicts (PyMOL error: loading mmCIF into existing object not supported)
cmd.load("E:\OneDrive - u-bordeaux.fr\Projects\Dopamine\images\multi.min.pdb")

# Calculate RMS stats for all atoms (+/- hydrogens)
apt_and_dop = calc_rms(
    "multi.min",
    selection="all",
    quiet=False,
    display_level=2,
    export_report=True,
    report_path="./.personal_use/reports"
)

apt_only = calc_rms(
    "multi.min",
    selection="not resn LDP",
    quiet=False,
    display_level=2,    
    export_report=True,
    report_path="./.personal_use/reports"
)

dop_only = calc_rms(   
    "multi.min",
    selection="resn LDP",
    quiet=False,
    display_level=2,    
    export_report=True,
    report_path="./.personal_use/reports"
)

apt_backbone = calc_rms(
    "multi.min",
    selection="backbone and not resn LDP",
    quiet=False,
    display_level=2,
    export_report=True,
    report_path="./.personal_use/reports"
)
