"""
Pairwise RMS deviation calculator for PyMOL.

Enhanced: calculate RMS values twice:
    - once with all selected atoms
    - once with hydrogens excluded

Usage from PyMOL:
        calculate_rms_stats obj_name, selection="all", quiet_py=False, verbose=True, display_level=2

display_level:
    0 -> show selection + average values table only
    1 -> show selection + average values table + average RMS per state
    2 -> show selection + average values table + average RMS per state + detailed per-state pair RMS

Returns:
        dict with keys:
            - "modes": dict with keys "all" and "no_hydrogen", each containing:
                    - "per_state_means": dict[state_index -> mean_rms or nan]
                    - "per_state_pairs": dict[state_index -> list of individual pair RMS values]
                    - "all_rms": list of individual RMS values
                    - "overall_mean": mean of all_rms (nan if empty)
                    - "overall_sd": stdev of all_rms (nan if <2 values)
            - "selection": resolved selection string used (base)
            - "atom_counts": dict for atom counts per mode (or None if unavailable)
            - "states": number of states inspected
            - "report_file": path to exported report if export_report True else None
            - "error": optional error message if something failed

If called with no arguments, the function returns a short help/usage message in a dict under the "help" key
and (unless quiet_py=True) prints it to the console.
"""
from typing import Dict, List, Any, Tuple, Optional
import math
import statistics
from pymol import cmd
from pathlib import Path
import datetime
import os


def calculate_rms_stats(
    obj_to_fit: Optional[str] = None,
    selection: str = "all",
    quiet: bool = False,
    display_level: int = 1,
    export_report: bool = False,
    report_path: str = ".",
) -> Dict[str, Any]:
    """
    Calculate RMS statistics across states for a PyMOL object/selection,
    reporting results for (1) all atoms and (2) excluding hydrogens.

    Parameters:
        obj_to_fit: object name or selection to fit. If omitted (None), the function
            returns a usage/help message instead of running.
        selection: additional selection string (default "all")
        quiet: when True suppress printing to console (default False)
        display_level: integer controlling printed output:
        0 -> only selection + average values table
        1 -> + average RMS per state
        2 -> + detailed per-state individual pair RMS report (default)
        export_report: when True export the report to a text file (same formatting as console)
        report_path: destination path; if a directory is provided the report will be saved there
             with a timestamped filename. Default is "." (current working directory).
    """
    # If called without specifying the required object, return usage/help info
    if obj_to_fit is None:
        help_lines = [
            "",
            "calculate_rms_stats usage:",
            "  calculate_rms_stats(obj_to_fit, selection='all', quiet=False, display_level=1, export_report=False, report_path='.')",
            "",
            "Arguments:",
            "  obj_to_fit   : (required) object name or selection to analyze. If omitted, this help is returned.",
            "  selection    : additional PyMOL selection (default: 'all')",
            "  quiet        : suppress printing when True",
            "  display_level: 0 -> summary only, 1 -> include per-state means, 2 -> full per-pair detail (default 1)",
            "  export_report: write report to file when True (default False)",
            "  report_path  : file path or directory for exported report, defaults to current working directory (default '.')",
            "",
            "Behavior:",
            "  The function computes RMS statistics twice: including all (selected) atoms, and excluding hydrogens.",
            "  Returns a dict with the analysis when run; when called with no obj_to_fit, returns this help:",
            ""
        ]
        help_text = "\n".join(help_lines)
        if not quiet:
            print(help_text)
        return {"help": help_text}

    if display_level not in (0, 1, 2):
        raise ValueError("display_level must be 0, 1, or 2")

    result: Dict[str, Any] = {
        "modes": {
            "all": {
                "per_state_means": {},
                "per_state_pairs": {},
                "all_rms": [],
                "overall_mean": math.nan,
                "overall_sd": math.nan,
            },
            "no_hydrogen": {
                "per_state_means": {},
                "per_state_pairs": {},
                "all_rms": [],
                "overall_mean": math.nan,
                "overall_sd": math.nan,
            },
        },
        "selection": None,
        "atom_counts": {"all": None, "no_hydrogen": None},
        "states": 0,
        "report_file": None,
        "error": None,
    }

    base_selection = f"({obj_to_fit}) and ({selection})"
    result["selection"] = base_selection

    def _p(*args, **kwargs):
        if not quiet:
            print(*args, **kwargs)

    try:
        states = cmd.count_states(obj_to_fit)
    except Exception as exc:
        err = f"Failed to determine states for '{obj_to_fit}': {exc}"
        result["error"] = err
        _p(err)
        return result

    result["states"] = states

    # Prepare both selection strings
    sel_all = base_selection
    sel_no_h = f"{base_selection} and (not hydro)"

    # Try to get atom counts for both selections (non-fatal)
    try:
        result["atom_counts"]["all"] = cmd.count_atoms(sel_all)
    except Exception as exc:
        _p(f"Warning: failed to count atoms for selection '{sel_all}': {exc}")
        result["atom_counts"]["all"] = None

    try:
        result["atom_counts"]["no_hydrogen"] = cmd.count_atoms(sel_no_h)
    except Exception as exc:
        _p(f"Warning: failed to count atoms for selection '{sel_no_h}': {exc}")
        result["atom_counts"]["no_hydrogen"] = None

    def compute_mode(selection_string: str) -> Tuple[Dict[int, float], List[float], Dict[int, List[float]]]:
        per_state_means: Dict[int, float] = {}
        all_rms: List[float] = []
        per_state_pairs: Dict[int, List[float]] = {}

        for state_idx in range(1, states + 1):
            try:
                rms_values = cmd.intra_fit(selection_string, state_idx, quiet=1)
            except Exception as exc:
                _p(f"Warning: intra_fit failed for state {state_idx} on selection '{selection_string}': {exc}")
                per_state_means[state_idx] = math.nan
                per_state_pairs[state_idx] = []
                continue

            if rms_values is None:
                rms_values = []
            elif isinstance(rms_values, (int, float)):
                rms_values = [float(rms_values)]
            else:
                try:
                    rms_values = [float(x) for x in rms_values]
                except Exception:
                    rms_values = []

            # Keep only non-negative (valid) values
            rms_filtered = [x for x in rms_values if x >= 0]

            per_state_pairs[state_idx] = rms_filtered[:]  # store list for this state

            if rms_filtered:
                m = statistics.mean(rms_filtered)
                per_state_means[state_idx] = m
                all_rms.extend(rms_filtered)
            else:
                per_state_means[state_idx] = math.nan

        return per_state_means, all_rms, per_state_pairs

    # Compute for both modes
    per_all, all_rms_all, per_pairs_all = compute_mode(sel_all)
    per_no_h, all_rms_no_h, per_pairs_no_h = compute_mode(sel_no_h)

    # Summarize and store
    def summarize(mode_dict: Dict[str, Any], per_state: Dict[int, float], rms_list: List[float],
              per_state_pairs: Dict[int, List[float]]):
        mode_dict["per_state_means"] = per_state
        mode_dict["per_state_pairs"] = per_state_pairs
        mode_dict["all_rms"] = rms_list
        if rms_list:
            mode_dict["overall_mean"] = statistics.mean(rms_list)
            mode_dict["overall_sd"] = statistics.stdev(rms_list) if len(rms_list) > 1 else math.nan
        else:
            mode_dict["overall_mean"] = math.nan
            mode_dict["overall_sd"] = math.nan

    summarize(result["modes"]["all"], per_all, all_rms_all, per_pairs_all)
    summarize(result["modes"]["no_hydrogen"], per_no_h, all_rms_no_h, per_pairs_no_h)

    # Build formatted report lines (returned as list of strings)
    def format_report_lines() -> List[str]:
        def fmt_num(x):
            return f"{x:.3f}" if (x is not None and not (isinstance(x, float) and math.isnan(x))) else "nan"

        lines: List[str] = []

        lines.append("")
        lines.append("*****************************")
        lines.append("* RMS Statistics Calculator *")
        lines.append("*****************************")
        lines.append("")

        if display_level == 0:
            display_level_str: str = "0 (average values only)"
        elif display_level == 1:
            display_level_str = "1 (global average and per-state means included)"
        else:
            display_level_str = "2 (full detail)"
        lines.append(f"The level of detail in this report has been set to => {display_level_str}")
        lines.append("")
        lines.append("############ Data selection ############")

        try:
            selection_provided = (selection is not None) and (str(selection).strip() != "") and (str(selection).strip().lower() != "all")
        except Exception:
            selection_provided = True  # conservative default

        if not selection_provided:
            lines.append("Note: no explicit selection was provided, so the script selected all atoms in the object.")
            lines.append("If you want a subset (for example only chain A), pass a selection string via the 'selection' argument, e.g.: selection='chain A'")
            lines.append("It is not necessary to explicitely exlude hydrogen atoms. RMS statistics are calculated both including and excluding hydrogens automatically.")
            lines.append("")
        lines.append("")
        lines.append("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        lines.append(f"Resolved selection: {base_selection}")
        natoms_all = result["atom_counts"]["all"]
        natoms_noh = result["atom_counts"]["no_hydrogen"]
        natoms_all_str = str(natoms_all) if natoms_all is not None else "?"
        natoms_noh_str = str(natoms_noh) if natoms_noh is not None else "?"
        lines.append(f"States inspected: {states}")
        lines.append("Atoms selected:")
        lines.append(f"  All: {natoms_all_str}")
        lines.append(f"  Heavy (no H): {natoms_noh_str}")
        lines.append("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        lines.append("")
        lines.append("########################################")
        lines.append("")

        # 2) Average value table (always shown)
        lines.append("======================== RMS Summary ========================")
        lines.append(f"{'Mode':<16} {'#atoms':>7} {'#values':>9} {'mean':>10} {'sd':>10}")
        lines.append("-" * 60)

        rows = [
            ("All atoms", result["atom_counts"]["all"], len(result["modes"]["all"]["all_rms"]),
             result["modes"]["all"]["overall_mean"], result["modes"]["all"]["overall_sd"]),
            ("Heavy atoms", result["atom_counts"]["no_hydrogen"], len(result["modes"]["no_hydrogen"]["all_rms"]),
             result["modes"]["no_hydrogen"]["overall_mean"], result["modes"]["no_hydrogen"]["overall_sd"]),
        ]

        for name, natoms, nvals, meanv, sdv in rows:
            natoms_str = str(natoms) if natoms is not None else "?"
            lines.append(f"{name:<16} {natoms_str:>7} {nvals:>9} {fmt_num(meanv):>10} {fmt_num(sdv):>10}")

        lines.append("-" * 60)

        # 3) Average RMS per state (display_level >= 1)
        if display_level >= 1:
            lines.append("Average RMS per state (showing NaN where unavailable):")
            lines.append("")
            lines.append("All atoms:")
            lines.append(f"{'State':>5} {'Mean RMS':>12} {'#pairs':>10}")
            lines.append("-" * 34)
            for state_idx in range(1, states + 1):
                m_all = result["modes"]["all"]["per_state_means"].get(state_idx, math.nan)
                n_pairs_all = len(result["modes"]["all"]["per_state_pairs"].get(state_idx, []))
                lines.append(f"{state_idx:5d} {fmt_num(m_all):>12} {n_pairs_all:>10d}")
            lines.append("")

            lines.append("Heavy atoms (no H):")
            lines.append(f"{'State':>5} {'Mean RMS':>12} {'#pairs':>10}")
            lines.append("-" * 34)
            for state_idx in range(1, states + 1):
                m_noh = result["modes"]["no_hydrogen"]["per_state_means"].get(state_idx, math.nan)
                n_pairs_noh = len(result["modes"]["no_hydrogen"]["per_state_pairs"].get(state_idx, []))
                lines.append(f"{state_idx:5d} {fmt_num(m_noh):>12} {n_pairs_noh:>10d}")
            lines.append("-" * 60)

        # 4) Detailed per-state individual pair RMS values (display_level == 2)
        if display_level >= 2:
            lines.append("Detailed per-state individual pair RMS values:")
            for state_idx in range(1, states + 1):
                lines.append("")
                lines.append(f"State {state_idx} - All atoms (mean = {fmt_num(result['modes']['all']['per_state_means'].get(state_idx, math.nan))})")
                lines.append(f"{'Pair #':>7} {'RMS':>12}")
                lines.append("-" * 22)
                pairs_all = result["modes"]["all"]["per_state_pairs"].get(state_idx, [])
                if pairs_all:
                    for i, val in enumerate(pairs_all, start=1):
                        lines.append(f"{i:7d} {fmt_num(val):>12}")
                else:
                    lines.append("   (no pair RMS values)")

                lines.append("")  # spacer
                lines.append(f"State {state_idx} - Heavy atoms (mean = {fmt_num(result['modes']['no_hydrogen']['per_state_means'].get(state_idx, math.nan))})")
                lines.append(f"{'Pair #':>7} {'RMS':>12}")
                lines.append("-" * 22)
                pairs_noh = result["modes"]["no_hydrogen"]["per_state_pairs"].get(state_idx, [])
                if pairs_noh:
                    for i, val in enumerate(pairs_noh, start=1):
                        lines.append(f"{i:7d} {fmt_num(val):>12}")
                else:
                    lines.append("   (no pair RMS values)")
                lines.append("-" * 60)
        else:
            if display_level == 0:
                lines.append("(Average RMS per state and detailed per-pair output suppressed; set display_level>=1 to enable per-state means or display_level==2 for full detail)")
            elif display_level == 1:
                lines.append("(Detailed per-pair output suppressed; set display_level==2 to enable)")

        lines.append("========================================================")
        lines.append("")

        return lines

    # Present to console (subject to quiet) and optionally export to file
    report_lines = format_report_lines()

    # Print to console using existing _p wrapper
    for ln in report_lines:
        _p(ln)

    # Export to file if requested (write even if quiet True)
    if export_report:
        try:
            target = Path(report_path)
            if target.exists() and target.is_dir():
                # sanitize object name and selection for a safe filename
                safe_obj = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in str(obj_to_fit))
                safe_obj = safe_obj[:30]  # keep reasonably short

                sel_str = "" if selection is None else str(selection)
                include_sel = sel_str.strip() != "" and sel_str.strip().lower() != "all"
                if include_sel:
                    safe_sel = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in sel_str)
                    # collapse repeated underscores and limit length
                    while "__" in safe_sel:
                        safe_sel = safe_sel.replace("__", "_")
                    safe_sel = safe_sel.strip("_")[:40]
                    sel_part = f"_sel_{safe_sel}"
                else:
                    sel_part = ""

                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                filename = f"rms_report_{safe_obj}{sel_part}_{timestamp}.txt"
                file_path = target / filename
            else:
                # treat provided path as file path (parent may not exist)
                file_path = target

            # Ensure parent directory exists
            if not file_path.parent.exists():
                file_path.parent.mkdir(parents=True, exist_ok=True)

            file_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
            result["report_file"] = str(file_path.resolve())
            if not quiet:
                _p(f"Report exported to: {result['report_file']}")
        except Exception as exc:
            err = f"Failed to write report to '{report_path}': {exc}"
            result["error"] = err if result.get("error") is None else f"{result['error']}; {err}"
            _p(err)
    else:
        _p("To export the report to a file, set export_report=True")
        _p("Set a specified path with the 'report_path' argument.")
        _p(f"In absence of a specified path, the report will be saved in the current working directory: {os.getcwd()}")

    return result


# Register the function as a PyMOL command
cmd.extend("calculate_rms_stats", calculate_rms_stats)
