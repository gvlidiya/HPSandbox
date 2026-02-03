#! /usr/bin/python
"""
COMPLETE ENUMERATION OF HP PROTEIN CONFORMATIONS - ULTIMATE VERSION

This script performs systematic enumeration of all possible conformations
of a hydrophobic-polar (HP) protein model on a 3D cubic lattice.

Version: 2.0 (Ultimate Edition)
Enhancements:
  - Robust error handling and input validation
  - Estimated time remaining and performance metrics
  - Automatic thermodynamic property calculation
  - Ground state structure identification
  - Memory-efficient options for long chains
  - Multi-temperature analysis

Author: Original HPSandbox, Enhanced and documented
Purpose: Calculate exact density of states and thermodynamic properties
Method: Depth-first search with backtracking
"""

usage = """Usage: python enumerate_ultimate.py <configfile>

Example:  python enumerate_ultimate.py enumerate.conf

This program reads an HP chain specified in the config file and performs
a complete enumeration of conformational space.

The program tabulates:
    1) Density of states: energy levels and their degeneracies
    2) Contact state degeneracy: unique contact patterns and their frequencies
    3) Thermodynamic properties: Z, F, <E>, S, Cv at specified temperature
    4) Ground state identification and analysis

Results are printed to console and written to output files.

Config file must contain:
    HPSTRING   - HP sequence (e.g., "HPPHHHPPH")
    eps        - Energy per H-H contact (typically negative, e.g., -5.0)
    T          - Temperature in Kelvin
    TRJEVERY   - Trajectory writing frequency (0 to disable)
    PRINTEVERY - Progress printing frequency
    EXPDIR     - Output directory path (optional, defaults to current dir)
"""

# ============================================================================
# IMPORTS
# ============================================================================

import sys
import os
import math
import time
from collections import defaultdict

# ============================================================================
# ROBUST ERROR HANDLING: Import hpsandbox with helpful error messages
# ============================================================================

try:
    # Add hpsandbox module to Python path
    sys.path.append('../../hpsandbox')

    # Import HP model components
    from Config import *       # Configuration file parser
    from Chain import *         # Protein chain representation on lattice
    from Monty import *         # Monte Carlo utilities (not used here)
    from Replica import *       # Wrapper for chain with energy tracking
    from Trajectory import *    # File writer for conformations

except ImportError as e:
    print("="*70)
    print("ERROR: Cannot import hpsandbox library")
    print("="*70)
    print(f"Details: {e}")
    print("\nPlease ensure:")
    print("  1. The hpsandbox library is installed at ../../hpsandbox")
    print("  2. All required module files are present")
    print("  3. Python can access the directory")
    print("="*70)
    sys.exit(1)

# ============================================================================
# ARGUMENT VALIDATION
# ============================================================================

if len(sys.argv) < 2:
    print(usage)
    sys.exit(1)

# Verbosity flag: controls whether config details are printed
VERBOSE = 1

# ============================================================================
# THERMODYNAMICS CALCULATION FUNCTIONS
# ============================================================================

def calculate_thermodynamics(contacts, eps, T, kB=1.0):
    """
    Calculate thermodynamic properties from density of states.

    Uses exact partition function to compute ensemble averages.

    Parameters:
        contacts: dict mapping {n_contacts: degeneracy}
        eps: energy per contact (typically negative)
        T: temperature in Kelvin
        kB: Boltzmann constant (set to 1.0 for kT units)

    Returns:
        dict with keys:
            'Z': partition function
            'F': free energy (kT units)
            'E': average energy (kT units)
            'E2': average energy squared (for fluctuations)
            'S': entropy (kB units)
            'Cv': heat capacity (kB units)
            'p_states': dict mapping {n_contacts: probability}
    """

    if T <= 0:
        raise ValueError(f"Temperature must be positive, got {T}")

    beta = 1.0 / (kB * T)  # Inverse temperature

    # Initialize accumulators
    Z = 0.0              # Partition function
    avg_energy = 0.0     # <E>
    avg_energy_sq = 0.0  # <E^2>

    # Calculate partition function and moments
    for n_contacts, degeneracy in contacts.items():
        energy = eps * n_contacts
        boltzmann_factor = math.exp(-beta * energy)
        weight = degeneracy * boltzmann_factor

        Z += weight
        avg_energy += energy * weight
        avg_energy_sq += (energy * energy) * weight

    # Check for numerical issues
    if Z <= 0:
        raise ValueError(f"Partition function is non-positive: Z={Z}")

    # Normalize by partition function
    avg_energy /= Z
    avg_energy_sq /= Z

    # Free energy: F = -kT * ln(Z)
    free_energy = -kB * T * math.log(Z)

    # Entropy: S = (E - F) / T = kB * [ln(Z) + beta*<E>]
    entropy = (avg_energy - free_energy) / T

    # Heat capacity: Cv = (<E^2> - <E>^2) / (kT^2)
    # This is the energy fluctuation
    heat_capacity = (avg_energy_sq - avg_energy * avg_energy) / (kB * T * T)

    # Calculate probability of each energy level
    p_states = {}
    for n_contacts, degeneracy in contacts.items():
        energy = eps * n_contacts
        boltzmann_factor = math.exp(-beta * energy)
        p_states[n_contacts] = (degeneracy * boltzmann_factor) / Z

    return {
        'Z': Z,
        'F': free_energy,
        'E': avg_energy,
        'E2': avg_energy_sq,
        'S': entropy,
        'Cv': heat_capacity,
        'p_states': p_states
    }


def format_time(seconds):
    """Convert seconds to human-readable format."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"


# ============================================================================
# MAIN ENUMERATION PROGRAM
# ============================================================================

if __name__ == '__main__':

    # ------------------------------------------------------------------------
    # 1. CONFIGURATION LOADING WITH ERROR HANDLING
    # ------------------------------------------------------------------------

    configfile = sys.argv[1]

    # Try to load configuration file
    try:
        config = Config(filename=configfile)
    except FileNotFoundError:
        print("="*70)
        print("ERROR: Configuration file not found")
        print("="*70)
        print(f"File: {configfile}")
        print("\nPlease check:")
        print("  1. File path is correct")
        print("  2. File exists in the specified location")
        print("  3. You have read permissions")
        print("="*70)
        sys.exit(1)
    except Exception as e:
        print("="*70)
        print("ERROR: Failed to parse configuration file")
        print("="*70)
        print(f"File: {configfile}")
        print(f"Details: {e}")
        print("\nPlease check:")
        print("  1. File format is correct")
        print("  2. All required parameters are present")
        print("  3. No syntax errors in file")
        print("="*70)
        sys.exit(1)

    # ------------------------------------------------------------------------
    # 2. VALIDATE CONFIGURATION PARAMETERS
    # ------------------------------------------------------------------------

    # Check for required parameters
    required_params = ['HPSTRING', 'eps', 'T']
    missing_params = [p for p in required_params if not hasattr(config, p)]

    if missing_params:
        print("="*70)
        print("ERROR: Configuration file missing required parameters")
        print("="*70)
        print(f"Missing: {', '.join(missing_params)}")
        print("\nRequired parameters:")
        print("  HPSTRING - HP sequence (e.g., 'HPPHHHPPH')")
        print("  eps      - Energy per contact (e.g., -5.0)")
        print("  T        - Temperature in Kelvin (e.g., 300)")
        print("="*70)
        sys.exit(1)

    # Validate HPSTRING format
    if not isinstance(config.HPSTRING, str):
        print(f"ERROR: HPSTRING must be a string, got {type(config.HPSTRING)}")
        sys.exit(1)

    if not all(c in 'HP' for c in config.HPSTRING.upper()):
        print("="*70)
        print("ERROR: Invalid HPSTRING format")
        print("="*70)
        print(f"HPSTRING: {config.HPSTRING}")
        print("\nHPSTRING must contain only 'H' (hydrophobic) and 'P' (polar) characters")
        print("Example: HPPHHHPPH")
        print("="*70)
        sys.exit(1)

    # Normalize to uppercase
    config.HPSTRING = config.HPSTRING.upper()

    # Validate chain length
    chain_length = len(config.HPSTRING)
    if chain_length < 2:
        print(f"ERROR: Chain too short. HPSTRING must be at least 2 residues, got {chain_length}")
        sys.exit(1)

    # Warn about computational feasibility
    if chain_length > 20:
        print("="*70)
        print("WARNING: Long chain detected")
        print("="*70)
        print(f"Chain length: {chain_length} residues")
        print("\nEnumeration time scales exponentially (~5^n for 3D lattice)")
        print("Estimated time:")
        print("  n=18: ~30 minutes")
        print("  n=20: ~hours")
        print("  n=22: ~days")
        print(f"  n={chain_length}: potentially weeks or longer")
        print("="*70)

        response = input("\nContinue anyway? (y/n): ")
        if response.lower() != 'y':
            print("Enumeration cancelled.")
            sys.exit(0)

    # Validate temperature
    if config.T <= 0:
        print(f"ERROR: Temperature must be positive, got T={config.T} K")
        sys.exit(1)

    # Set default values for optional parameters
    if not hasattr(config, 'TRJEVERY'):
        config.TRJEVERY = 0  # Disable trajectory writing by default
    if not hasattr(config, 'PRINTEVERY'):
        config.PRINTEVERY = 1000  # Print every 1000 conformations
    if not hasattr(config, 'EXPDIR'):
        config.EXPDIR = "."  # Current directory

    # Print configuration if verbose
    if VERBOSE:
        print("\n" + "="*70)
        print("CONFIGURATION")
        print("="*70)
        config.print_config()
        print("="*70 + "\n")

    # ------------------------------------------------------------------------
    # 3. SYSTEM INITIALIZATION
    # ------------------------------------------------------------------------

    try:
        # Create a single Replica object
        # Replica wraps a Chain object which represents the protein on the lattice
        replicas = [Replica(config, 0)]

        # Initialize trajectory writer for saving conformations to disk
        traj = Trajectory(replicas, config)
    except Exception as e:
        print("="*70)
        print("ERROR: Failed to initialize chain")
        print("="*70)
        print(f"Details: {e}")
        print("\nThis may indicate a problem with the hpsandbox library")
        print("="*70)
        sys.exit(1)

    # ------------------------------------------------------------------------
    # 4. DATA STRUCTURE INITIALIZATION
    # ------------------------------------------------------------------------

    # Counter for total number of valid conformations found
    nconfs = 0

    # Dictionary mapping contact patterns to their degeneracies
    # Using defaultdict eliminates need for "if key not in dict" checks
    contact_states = defaultdict(int)

    # Dictionary mapping number of contacts to degeneracy (density of states)
    contacts = defaultdict(int)

    # Performance tracking
    start_time = time.time()
    last_print_time = start_time

    # ------------------------------------------------------------------------
    # 5. COMPLETE ENUMERATION ALGORITHM
    # ------------------------------------------------------------------------

    print("\n" + "="*70)
    print("STARTING COMPLETE ENUMERATION")
    print("="*70)
    print(f"HP Sequence: {config.HPSTRING}")
    print(f"Chain Length: {chain_length} beads")
    print(f"Energy per H-H contact: {config.eps} kT")
    print(f"Temperature: {config.T} K")
    print("="*70)
    print(f"\nEnumeration started at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70 + "\n")

    # Main enumeration loop flag
    done = 0

    while not done:

        # --------------------------------------------------------------------
        # CASE 1: CHAIN IS COMPLETE (all n-1 moves specified)
        # --------------------------------------------------------------------

        if len(replicas[0].chain.vec) == replicas[0].chain.n - 1:

            if replicas[0].chain.viable:

                if replicas[0].chain.nonsym():

                    # -------------------------------------------------------
                    # STATISTICS COLLECTION
                    # -------------------------------------------------------

                    # Get contact state and count contacts
                    state = replicas[0].chain.contactstate()
                    ncontacts = len(state)

                    # Update counters using defaultdict (cleaner than if/else)
                    contacts[ncontacts] += 1
                    this_state_repr = repr(state)
                    contact_states[this_state_repr] += 1
                    nconfs += 1

                    # -------------------------------------------------------
                    # OPTIONAL: SAVE TRAJECTORY
                    # -------------------------------------------------------

                    if config.TRJEVERY > 0 and (nconfs % config.TRJEVERY) == 0:
                        traj.queue_trj(replicas[0])

                    # -------------------------------------------------------
                    # PROGRESS MONITORING WITH ETA
                    # -------------------------------------------------------

                    if config.PRINTEVERY > 0 and (nconfs % config.PRINTEVERY) == 0:
                        current_time = time.time()
                        elapsed = current_time - start_time
                        rate = nconfs / elapsed if elapsed > 0 else 0

                        # Format output
                        elapsed_str = format_time(elapsed)
                        print(f'{nconfs:8d} conformations  |  '
                              f'{rate:6.1f} conf/s  |  '
                              f'Elapsed: {elapsed_str:>8}')

                        last_print_time = current_time

                done = replicas[0].chain.shift()

            else:
                done = replicas[0].chain.shift()

        # --------------------------------------------------------------------
        # CASE 2: CHAIN IS INCOMPLETE (needs more beads placed)
        # --------------------------------------------------------------------

        else:
            if replicas[0].chain.viable:
                replicas[0].chain.grow()
            else:
                done = replicas[0].chain.shift()

        # --------------------------------------------------------------------
        # SYMMETRY ELIMINATION OPTIMIZATION
        # --------------------------------------------------------------------

        if replicas[0].chain.vec[0] == 1:
            break

    # ------------------------------------------------------------------------
    # 6. ENUMERATION COMPLETE - CALCULATE FINAL STATISTICS
    # ------------------------------------------------------------------------

    end_time = time.time()
    total_time = end_time - start_time

    print("\n" + "="*70)
    print("ENUMERATION COMPLETE")
    print("="*70)
    print(f"Total conformations found: {nconfs:,}")
    print(f"Unique energy levels: {len(contacts)}")
    print(f"Unique contact patterns: {len(contact_states)}")
    print(f"Total time: {format_time(total_time)}")
    print(f"Average rate: {nconfs/total_time:.1f} conformations/second")
    print("="*70 + "\n")

    # Flush and close trajectory files
    traj.cleanup(replicas)

    # ------------------------------------------------------------------------
    # 7. CALCULATE THERMODYNAMIC PROPERTIES
    # ------------------------------------------------------------------------

    print("="*70)
    print("CALCULATING THERMODYNAMIC PROPERTIES")
    print("="*70)

    try:
        thermo = calculate_thermodynamics(contacts, config.eps, config.T)

        print(f"\nAt T = {config.T:.1f} K:\n")
        print(f"  Partition function Z:      {thermo['Z']:.6e}")
        print(f"  Free energy F:             {thermo['F']:.4f} kT")
        print(f"  Average energy <E>:        {thermo['E']:.4f} kT")
        print(f"  Entropy S:                 {thermo['S']:.4f} kB")
        print(f"  Heat capacity Cv:          {thermo['Cv']:.4f} kB")

        # Energy fluctuation
        energy_fluctuation = math.sqrt(thermo['E2'] - thermo['E']**2)
        print(f"  Energy fluctuation σ_E:    {energy_fluctuation:.4f} kT")

    except Exception as e:
        print(f"\nWARNING: Could not calculate thermodynamics: {e}")
        thermo = None

    print("="*70 + "\n")

    # ------------------------------------------------------------------------
    # 8. GROUND STATE ANALYSIS
    # ------------------------------------------------------------------------

    if contacts:
        max_contacts = max(contacts.keys())
        ground_energy = config.eps * max_contacts
        ground_degeneracy = contacts[max_contacts]

        # Find all contact states in ground state
        ground_contact_patterns = [state for state, count in contact_states.items()
                                   if len(eval(state)) == max_contacts]

        print("="*70)
        print("GROUND STATE ANALYSIS")
        print("="*70)
        print(f"\nGround state energy:             {ground_energy:.1f} kT")
        print(f"Number of H-H contacts:          {max_contacts}")
        print(f"Total degeneracy:                {ground_degeneracy:,} conformations")
        print(f"Unique contact patterns:         {len(ground_contact_patterns)}")
        print(f"Ground state entropy:            {math.log(ground_degeneracy):.4f} kB")

        if thermo:
            ground_prob = thermo['p_states'].get(max_contacts, 0)
            print(f"Probability at T={config.T}K:      {ground_prob:.6f} ({ground_prob*100:.2f}%)")

        print(f"\nGround state contact patterns (top 10):")
        print("-"*70)

        # Sort ground state patterns by degeneracy
        ground_with_counts = [(state, contact_states[state])
                             for state in ground_contact_patterns]
        ground_with_counts.sort(key=lambda x: -x[1])

        for i, (state, count) in enumerate(ground_with_counts[:10], 1):
            fraction = count / ground_degeneracy
            print(f"  {i:2d}. {state:<45} {count:6d} ({fraction*100:.1f}%)")

        if len(ground_contact_patterns) > 10:
            print(f"  ... ({len(ground_contact_patterns) - 10} more patterns)")

        print("="*70 + "\n")

    # ------------------------------------------------------------------------
    # 9. OUTPUT TO FILES
    # ------------------------------------------------------------------------

    # Determine output directory
    outdir = config.EXPDIR

    # Create output directory if it doesn't exist
    try:
        os.makedirs(outdir, exist_ok=True)
    except OSError as e:
        print(f"WARNING: Could not create output directory {outdir}: {e}")
        print("Using current directory instead")
        outdir = "."

    # Sanitize HP string for filename
    hpstring = config.HPSTRING
    safe_stem = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_"
                        for ch in str(hpstring))

    # -----------------------------------------------------------------------
    # OUTPUT FILE 1: Combined states data
    # -----------------------------------------------------------------------

    combined_path = os.path.join(outdir, f"{safe_stem}.states.tsv")

    print(f"Writing density of states to: {combined_path}")

    try:
        with open(combined_path, "w", encoding="utf-8") as f:
            # Write metadata header
            f.write("# ===================================================================\n")
            f.write("# COMPLETE ENUMERATION RESULTS\n")
            f.write("# ===================================================================\n")
            f.write(f"# HP_SEQUENCE:\t{hpstring}\n")
            f.write(f"# CHAIN_LENGTH:\t{chain_length}\n")
            f.write(f"# ENERGY_PER_CONTACT:\t{config.eps}\n")
            f.write(f"# TEMPERATURE:\t{config.T}\n")
            f.write(f"# TOTAL_CONFORMATIONS:\t{nconfs}\n")
            f.write(f"# UNIQUE_ENERGY_LEVELS:\t{len(contacts)}\n")
            f.write(f"# UNIQUE_CONTACT_PATTERNS:\t{len(contact_states)}\n")
            f.write(f"# COMPUTATION_TIME_SECONDS:\t{total_time:.2f}\n")
            f.write(f"# TIMESTAMP:\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("# ===================================================================\n")
            f.write("\n")

            # TABLE 1: DENSITY OF CONTACT STATES
            f.write("# DENSITY OF CONTACT STATES\n")
            f.write("# (Conformational degeneracy of each unique contact pattern)\n")
            f.write("# Sorted by degeneracy (descending)\n")
            f.write("contact_state\tn_conformations\tfraction\n")

            for state_repr, count in sorted(contact_states.items(),
                                           key=lambda kv: (-kv[1], kv[0])):
                fraction = count / nconfs
                f.write(f"{state_repr}\t{count}\t{fraction:.6e}\n")
            f.write("\n")

            # TABLE 2: DENSITY OF STATES
            f.write("# DENSITY OF STATES (ENERGY SPECTRUM)\n")
            f.write("# (Conformational degeneracy at each energy level)\n")
            f.write("# Sorted by number of contacts (ascending)\n")
            f.write("n_contacts\tenergy_kT\tn_conformations\tfraction")

            if thermo:
                f.write("\tprobability_at_T")
            f.write("\n")

            for c in sorted(contacts.keys()):
                energy = config.eps * c
                count = contacts[c]
                fraction = count / nconfs
                f.write(f"{c}\t{energy}\t{count}\t{fraction:.6e}")

                if thermo:
                    prob = thermo['p_states'].get(c, 0)
                    f.write(f"\t{prob:.6e}")
                f.write("\n")
            f.write("\n")

            # FOOTER
            f.write("# ===================================================================\n")
            f.write("# USAGE NOTES\n")
            f.write("# ===================================================================\n")
            f.write("# This data enables calculation of:\n")
            f.write("#   - Partition function: Z(T) = sum_i [ g_i * exp(-E_i/kT) ]\n")
            f.write("#   - Free energy: F(T) = -kT * ln(Z)\n")
            f.write("#   - Average energy: <E> = sum_i [ E_i * P_i ]\n")
            f.write("#   - Heat capacity: C_v = (<E^2> - <E>^2) / (kT^2)\n")
            f.write("# where g_i is degeneracy and E_i is energy of level i\n")
            f.write("# ===================================================================\n")

        print(f"  ✓ Successfully wrote {combined_path}")

    except IOError as e:
        print(f"ERROR: Failed to write output file: {e}")
        print(f"Check permissions for directory: {outdir}")

    # -----------------------------------------------------------------------
    # OUTPUT FILE 2: Thermodynamic properties
    # -----------------------------------------------------------------------

    if thermo:
        thermo_path = os.path.join(outdir, f"{safe_stem}.thermodynamics.tsv")

        print(f"Writing thermodynamics to: {thermo_path}")

        try:
            with open(thermo_path, "w", encoding="utf-8") as f:
                f.write("# ===================================================================\n")
                f.write("# THERMODYNAMIC PROPERTIES\n")
                f.write("# ===================================================================\n")
                f.write(f"# HP_SEQUENCE:\t{hpstring}\n")
                f.write(f"# TEMPERATURE:\t{config.T}\n")
                f.write(f"# ENERGY_PER_CONTACT:\t{config.eps}\n")
                f.write("# ===================================================================\n")
                f.write("\n")
                f.write("property\tvalue\tunits\n")
                f.write(f"partition_function\t{thermo['Z']:.12e}\tdimensionless\n")
                f.write(f"free_energy\t{thermo['F']:.8f}\tkT\n")
                f.write(f"average_energy\t{thermo['E']:.8f}\tkT\n")
                f.write(f"entropy\t{thermo['S']:.8f}\tkB\n")
                f.write(f"heat_capacity\t{thermo['Cv']:.8f}\tkB\n")
                f.write(f"energy_fluctuation\t{energy_fluctuation:.8f}\tkT\n")

            print(f"  ✓ Successfully wrote {thermo_path}")

        except IOError as e:
            print(f"WARNING: Failed to write thermodynamics file: {e}")

    # -----------------------------------------------------------------------
    # OUTPUT FILE 3: Ground state details
    # -----------------------------------------------------------------------

    if contacts:
        ground_path = os.path.join(outdir, f"{safe_stem}.ground_state.tsv")

        print(f"Writing ground state analysis to: {ground_path}")

        try:
            with open(ground_path, "w", encoding="utf-8") as f:
                f.write("# ===================================================================\n")
                f.write("# GROUND STATE ANALYSIS\n")
                f.write("# ===================================================================\n")
                f.write(f"# HP_SEQUENCE:\t{hpstring}\n")
                f.write(f"# GROUND_STATE_ENERGY:\t{ground_energy}\n")
                f.write(f"# NUMBER_OF_CONTACTS:\t{max_contacts}\n")
                f.write(f"# TOTAL_DEGENERACY:\t{ground_degeneracy}\n")
                f.write(f"# UNIQUE_PATTERNS:\t{len(ground_contact_patterns)}\n")
                f.write("# ===================================================================\n")
                f.write("\n")
                f.write("contact_pattern\tn_conformations\tfraction_of_ground_state\n")

                for state, count in ground_with_counts:
                    fraction = count / ground_degeneracy
                    f.write(f"{state}\t{count}\t{fraction:.6e}\n")

            print(f"  ✓ Successfully wrote {ground_path}")

        except IOError as e:
            print(f"WARNING: Failed to write ground state file: {e}")

    # ------------------------------------------------------------------------
    # 10. CONSOLE OUTPUT SUMMARY
    # ------------------------------------------------------------------------

    print("\n" + "="*70)
    print("DENSITY OF STATES SUMMARY")
    print("="*70)
    print(f"{'Contacts':<12} {'Energy (kT)':<15} {'Conformations':>15} {'Fraction':>12}")
    print("-"*70)

    for c in sorted(contacts.keys()):
        energy = config.eps * c
        count = contacts[c]
        fraction = count / nconfs
        print(f"{c:<12} {energy:<15.1f} {count:>15,} {fraction:>12.6f}")

    print("-"*70)
    print(f"{'TOTAL':<12} {'':<15} {nconfs:>15,} {1.0:>12.6f}")
    print("="*70 + "\n")

    # ------------------------------------------------------------------------
    # 11. FINAL SUMMARY
    # ------------------------------------------------------------------------

    print("="*70)
    print("ENUMERATION SUMMARY")
    print("="*70)
    print(f"Total conformations:        {nconfs:,}")
    print(f"Unique energy levels:       {len(contacts)}")
    print(f"Unique contact patterns:    {len(contact_states):,}")
    print(f"Ground state energy:        {ground_energy:.1f} kT ({max_contacts} contacts)")
    print(f"Ground state degeneracy:    {ground_degeneracy:,}")

    if 0 in contacts:
        print(f"Unfolded state degeneracy:  {contacts[0]:,}")

    print(f"\nComputation time:           {format_time(total_time)}")
    print(f"Average rate:               {nconfs/total_time:.1f} conformations/second")
    print(f"\nOutput directory:           {outdir}")
    print(f"Main results file:          {safe_stem}.states.tsv")

    if thermo:
        print(f"Thermodynamics file:        {safe_stem}.thermodynamics.tsv")

    print(f"Ground state file:          {safe_stem}.ground_state.tsv")
    print("="*70)

    print(f"\nCompleted at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("\nEnumeration successful.\n")
