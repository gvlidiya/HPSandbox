#! /usr/bin/python
"""
COMPLETE ENUMERATION OF HP PROTEIN CONFORMATIONS ON 2D SQUARE LATTICE

This script performs systematic enumeration of all possible conformations
of a hydrophobic-polar (HP) protein model on a 2D square lattice.

Author: Original HPSandbox, Enhanced with detailed comments
Purpose: Calculate exact density of states and contact state degeneracies
Method: Depth-first search with backtracking on 2D square lattice
Lattice: 2D square lattice with 4 cardinal directions (Up, Right, Down, Left)
"""

usage = """Usage: python enumerate_2D.py <configfile>

Try:  enumerate_2D.py enumerate.conf

This program reads an HP chain specified in the config file and performs
a complete enumeration of conformational space on a 2D square lattice.

The program tabulates:
    1) Density of states: energy levels and their degeneracies
    2) Contact state degeneracy: unique contact patterns and their frequencies

Results are printed to console and written to output files.
"""

# ============================================================================
# IMPORTS
# ============================================================================

import sys
import os
from collections import defaultdict

# Add hpsandbox module to Python path
# This allows importing the HP protein simulation library
sys.path.append('../../')

# Import HP model components
from hpsandbox import Chain, Config, Monty, Replica, Trajectory
# ============================================================================
# ARGUMENT VALIDATION
# ============================================================================

# Check for config file argument before entering main block
if len(sys.argv) < 2:
    print(usage)
    sys.exit(1)

# Verbosity flag: controls whether config details are printed
VERBOSE = 1

# ============================================================================
# MAIN ENUMERATION PROGRAM
# ============================================================================

if __name__ == '__main__':

    # ------------------------------------------------------------------------
    # 1. SYSTEM INITIALIZATION
    # ------------------------------------------------------------------------

    # Re-check arguments (defensive programming)
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(1)

    # Load configuration file containing:
    #   - HPSTRING: sequence like "HPPHHHPPH" (H=hydrophobic, P=polar)
    #   - eps: energy per H-H contact (typically negative, e.g., -5.0)
    #   - T: temperature in Kelvin
    #   - TRJEVERY: trajectory writing frequency
    #   - PRINTEVERY: progress printing frequency
    #   - EXPDIR: output directory path
    configfile = sys.argv[1]
    config = Config(filename=configfile)

    # Print configuration to console if verbose mode enabled
    if VERBOSE:
        config.print_config()

    # Create a single Replica object
    # Replica wraps a Chain object which represents the protein on the 2D square lattice
    # The Chain starts with move vector [0, 0, 0, ..., 0] (all moves "up" initially)
    # Move encoding: 0=Up, 1=Right, 2=Down, 3=Left
    replicas = [Replica(config, 0)]

    # Initialize trajectory writer for saving conformations to disk
    # This allows visualization and post-analysis of individual structures
    traj = Trajectory(replicas, config)

    # ------------------------------------------------------------------------
    # 2. DATA STRUCTURE INITIALIZATION
    # ------------------------------------------------------------------------

    # Counter for total number of valid conformations found
    nconfs = 0

    # Dictionary mapping contact patterns to their degeneracies
    # Key: string representation of contact list, e.g., "[(0, 9), (3, 7)]"
    # Value: number of conformations with that exact contact pattern
    # Using defaultdict eliminates need for "if key not in dict" checks
    contact_states = defaultdict(int)

    # Dictionary mapping number of contacts to degeneracy (density of states)
    # Key: number of H-H contacts (integer)
    # Value: number of conformations with that many contacts
    contacts = defaultdict(int)

    # ------------------------------------------------------------------------
    # 3. COMPLETE ENUMERATION ALGORITHM (2D SQUARE LATTICE)
    # ------------------------------------------------------------------------
    #
    # This implements a depth-first search through conformational space:
    #   - GROW: Extend chain by placing next bead in one of 4 lattice directions
    #   - CHECK: Verify conformation is viable (no overlaps, connected)
    #   - COUNT: If complete and valid, record statistics
    #   - SHIFT: Backtrack and try alternative placements
    #
    # The algorithm systematically explores ALL possible ways to fold the
    # chain on the 2D square lattice, subject to constraints:
    #   - Self-avoidance: no two beads occupy same (x,y) site
    #   - Connectivity: consecutive beads must be adjacent (Manhattan distance = 1)
    #   - No duplicates: rotational/mirror symmetries eliminated
    #
    # 2D SQUARE LATTICE SPECIFICS:
    #   - 4 possible moves per step: Up(0), Right(1), Down(2), Left(3)
    #   - Coordinates are 2D tuples: (x, y)
    #   - First bead fixed at origin: (0, 0)
    #   - Conformational space scales as ~3^n (each step has ~3 choices on average)
    # ------------------------------------------------------------------------

    print("\n" + "="*70)
    print("STARTING COMPLETE ENUMERATION - 2D SQUARE LATTICE")
    print("="*70)
    print(f"HP Sequence: {config.HPSTRING}")
    print(f"Chain Length: {replicas[0].chain.n} beads")
    print(f"Energy per H-H contact: {config.eps} kT")
    print(f"Temperature: {config.T} K")
    print(f"Lattice Type: 2D square lattice")
    print(f"Move Set: 4 directions (Up=0, Right=1, Down=2, Left=3)")
    print("="*70 + "\n")

    # Main enumeration loop flag
    # done=0 means more conformations to explore, done=1 means exhausted
    done = 0

    while not done:

        # --------------------------------------------------------------------
        # CASE 1: CHAIN IS COMPLETE (all n-1 moves specified)
        # --------------------------------------------------------------------
        # For a chain of n beads, we need n-1 lattice moves
        # When vec has length n-1, the chain is fully specified

        if len(replicas[0].chain.vec) == replicas[0].chain.n - 1:

            # Check if this complete conformation is viable
            # viable = True if:
            #   - No two beads occupy the same (x,y) lattice site (self-avoidance)
            #   - All beads form a connected path (no breaks)
            #   - All bond lengths equal 1 lattice unit (implicit in move encoding)

            if replicas[0].chain.viable:

                # Check if this conformation is non-symmetric
                # nonsym() returns True if this is NOT a rotation or reflection
                # of a previously counted conformation
                # Symmetry rule for 2D square lattice:
                #   - First move must be Up (0)
                #   - First turn must be Right (1) or continue Up (0)
                # This ensures we count each unique structure only once

                if replicas[0].chain.nonsym():

                    # -------------------------------------------------------
                    # STATISTICS COLLECTION
                    # -------------------------------------------------------

                    # Get the contact state: list of (i, j) pairs where:
                    #   - Both beads i and j are hydrophobic (H)
                    #   - |i - j| >= 2 (not bonded neighbors along chain)
                    #   - Beads are adjacent on 2D lattice (Manhattan distance = 1)
                    #     i.e., |x_i - x_j| + |y_i - y_j| = 1
                    # Example: [(0, 9), (3, 7)] means bead 0 touches bead 9,
                    #          and bead 3 touches bead 7 on the lattice
                    state = replicas[0].chain.contactstate()

                    # Count total number of H-H contacts
                    # This determines the energy: E = eps × ncontacts
                    ncontacts = len(state)

                    # Increment degeneracy counter for this energy level
                    # Using defaultdict, this automatically initializes to 0
                    contacts[ncontacts] += 1

                    # Convert contact state to string for use as dictionary key
                    # repr() gives a string like "[(0, 9), (3, 7)]"
                    this_state_repr = repr(state)

                    # Increment degeneracy counter for this specific contact pattern
                    contact_states[this_state_repr] += 1

                    # Increment total conformation counter
                    nconfs += 1

                    # -------------------------------------------------------
                    # OPTIONAL: SAVE TRAJECTORY
                    # -------------------------------------------------------
                    # Periodically save conformations to file for later visualization
                    # config.TRJEVERY controls frequency (e.g., save every 1000th)
                    if config.TRJEVERY > 0 and (nconfs % config.TRJEVERY) == 0:
                        traj.queue_trj(replicas[0])

                    # -------------------------------------------------------
                    # PROGRESS MONITORING
                    # -------------------------------------------------------
                    # Print progress to console so user knows it's working
                    # Shows conformation count and current move vector
                    if config.PRINTEVERY > 0 and (nconfs % config.PRINTEVERY) == 0:
                        print(f'{nconfs:6d} conformations  |  Current vector: {replicas[0].chain.vec}')

                # After counting (or skipping if symmetric), shift to next conformation
                # shift() retracts the last move and tries the next alternative
                # In 2D: cycles through Up(0) → Right(1) → Down(2) → Left(3)
                # Returns True when this branch is exhausted
                done = replicas[0].chain.shift()

            else:
                # Chain is complete but not viable (has overlaps)
                # Skip this conformation and shift to next alternative
                done = replicas[0].chain.shift()

        # --------------------------------------------------------------------
        # CASE 2: CHAIN IS INCOMPLETE (needs more beads placed)
        # --------------------------------------------------------------------
        else:
            # Check if current partial chain is viable
            # This pruning strategy saves time: if the partial chain already
            # has overlaps, no point extending it further

            if replicas[0].chain.viable:
                # Chain so far is valid, so grow it by one more bead
                # grow() adds next bead in direction 0 (Up) initially
                # The shift() method will cycle through all 4 directions
                replicas[0].chain.grow()
            else:
                # Partial chain already invalid (has overlap)
                # No point growing further, so backtrack immediately
                done = replicas[0].chain.shift()

        # --------------------------------------------------------------------
        # SYMMETRY ELIMINATION OPTIMIZATION (2D SPECIFIC)
        # --------------------------------------------------------------------
        # The move vector's first element (vec[0]) encodes the first move
        # For 2D square lattice with 4-fold rotational symmetry:
        #   - vec[0] = 0 (Up): explore this hemisphere
        #   - vec[0] = 1 (Right): rotations of vec[0]=0 cases
        #   - vec[0] = 2 (Down): rotations of vec[0]=0 cases
        #   - vec[0] = 3 (Left): rotations of vec[0]=0 cases
        # When vec[0] reaches 1, we've explored all unique conformations
        # due to the 4-fold rotational symmetry of the square lattice
        # This optimization reduces computation by ~4x

        if replicas[0].chain.vec[0] == 1:
            break

    # ------------------------------------------------------------------------
    # 4. FINALIZATION
    # ------------------------------------------------------------------------

    print("\n" + "="*70)
    print("ENUMERATION COMPLETE")
    print("="*70)
    print(f"Total conformations found: {nconfs}")
    print(f"Unique energy levels: {len(contacts)}")
    print(f"Unique contact patterns: {len(contact_states)}")
    print("="*70 + "\n")

    # Flush and close trajectory files
    # This ensures all buffered data is written to disk
    traj.cleanup(replicas)

    # ------------------------------------------------------------------------
    # 5. OUTPUT TO FILES
    # ------------------------------------------------------------------------

    # Determine output directory (from config or use current directory)
    try:
        outdir = config.EXPDIR
    except AttributeError:
        outdir = "."

    # Create output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)

    # Get HP string for filename (with fallback if not defined)
    try:
        hpstring = config.HPSTRING
    except AttributeError:
        hpstring = "UNKNOWN_HPSTRING"

    # Sanitize HP string to create safe filename
    # Replace any non-alphanumeric characters (except dash/underscore) with underscore
    # This prevents issues with spaces, slashes, etc. in filenames
    safe_stem = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_"
                        for ch in str(hpstring))

    # -----------------------------------------------------------------------
    # OUTPUT FILE: Combined density of states and contact state degeneracy
    # -----------------------------------------------------------------------

    combined_path = os.path.join(outdir, f"{safe_stem}.states_2D.tsv")

    print(f"Writing results to: {combined_path}\n")

    with open(combined_path, "w", encoding="utf-8") as f:
        # Write metadata header
        f.write("# ===================================================================\n")
        f.write("# COMPLETE ENUMERATION RESULTS - 2D SQUARE LATTICE\n")
        f.write("# ===================================================================\n")
        f.write(f"# HP_SEQUENCE:\t{hpstring}\n")
        f.write(f"# CHAIN_LENGTH:\t{replicas[0].chain.n}\n")
        f.write(f"# LATTICE_TYPE:\t2D_SQUARE\n")
        f.write(f"# LATTICE_DIMENSIONS:\t2\n")
        f.write(f"# MOVE_DIRECTIONS:\t4 (Up=0, Right=1, Down=2, Left=3)\n")
        f.write(f"# COORDINATE_SYSTEM:\t(x, y) tuples\n")
        f.write(f"# ENERGY_PER_CONTACT:\t{config.eps}\n")
        f.write(f"# TEMPERATURE:\t{config.T}\n")
        f.write(f"# TOTAL_CONFORMATIONS:\t{nconfs}\n")
        f.write(f"# UNIQUE_ENERGY_LEVELS:\t{len(contacts)}\n")
        f.write(f"# UNIQUE_CONTACT_PATTERNS:\t{len(contact_states)}\n")
        f.write("# ===================================================================\n")
        f.write("\n")

        # ---------------------------------------------------------------
        # TABLE 1: DENSITY OF CONTACT STATES
        # ---------------------------------------------------------------
        # Shows each unique contact pattern and how many conformations have it
        # Sorted by degeneracy (most common patterns first)

        f.write("# DENSITY OF CONTACT STATES\n")
        f.write("# (Conformational degeneracy of each unique contact pattern)\n")
        f.write("# Format: contact_state <TAB> n_conformations\n")
        f.write("# Sorted by degeneracy (descending)\n")
        f.write("contact_state\tn_conformations\n")

        # Sort by count (descending), then by contact state string (for ties)
        for state_repr, count in sorted(contact_states.items(),
                                       key=lambda kv: (-kv[1], kv[0])):
            f.write(f"{state_repr}\t{count}\n")
        f.write("\n")

        # ---------------------------------------------------------------
        # TABLE 2: DENSITY OF STATES (ENERGY SPECTRUM)
        # ---------------------------------------------------------------
        # Shows how many conformations exist at each energy level
        # This is the fundamental quantity for statistical mechanics

        f.write("# DENSITY OF STATES (ENERGY SPECTRUM)\n")
        f.write("# (Conformational degeneracy at each energy level)\n")
        f.write("# Format: n_contacts <TAB> energy_kT <TAB> n_conformations\n")
        f.write("# Sorted by number of contacts (ascending)\n")
        f.write("n_contacts\tenergy_kT\tn_conformations\n")

        # Sort by number of contacts (ascending = increasing energy)
        for c, count in sorted(contacts.items(), key=lambda kv: kv[0]):
            # Energy formula: E = eps × n_contacts
            # eps is typically negative (e.g., -5.0), so more contacts = lower energy
            energy = config.eps * c
            f.write(f"{c}\t{energy}\t{count}\n")
        f.write("\n")

        # ---------------------------------------------------------------
        # FOOTER: How to use this data
        # ---------------------------------------------------------------
        f.write("# ===================================================================\n")
        f.write("# USAGE NOTES\n")
        f.write("# ===================================================================\n")
        f.write("# This data enables calculation of:\n")
        f.write("#   - Partition function: Z(T) = sum_i [ g_i * exp(-E_i/kT) ]\n")
        f.write("#   - Free energy: F(T) = -kT * ln(Z)\n")
        f.write("#   - Average energy: <E> = sum_i [ E_i * P_i ]\n")
        f.write("#   - Heat capacity: C_v = d<E>/dT\n")
        f.write("#   - Probability of each state: P_i = g_i * exp(-E_i/kT) / Z\n")
        f.write("# where g_i is degeneracy and E_i is energy of state i\n")
        f.write("#\n")
        f.write("# 2D SQUARE LATTICE SPECIFIC NOTES:\n")
        f.write("#   - Contact detection: Manhattan distance = 1 in 2D plane\n")
        f.write("#   - Scaling: Conformational space ~ 3^n (approximately)\n")
        f.write("#   - Coordinates: (x, y) tuples representing positions on lattice\n")
        f.write("#   - Symmetry: 4-fold rotational + reflection symmetries eliminated\n")
        f.write("# ===================================================================\n")

    # ------------------------------------------------------------------------
    # 6. CONSOLE OUTPUT (for immediate inspection)
    # ------------------------------------------------------------------------

    print("="*70)
    print("DENSITY OF CONTACT STATES")
    print("="*70)
    print(f"{'Contact State':<50} {'Count':>15}")
    print("-"*70)

    # Show top 20 most common contact states (or all if fewer than 20)
    sorted_states = sorted(contact_states.items(), key=lambda kv: -kv[1])
    display_count = min(20, len(sorted_states))

    for i, (state, count) in enumerate(sorted_states[:display_count]):
        # Truncate very long contact state strings for display
        state_str = state if len(state) <= 45 else state[:42] + "..."
        print(f"{state_str:<50} {count:>15}")

    if len(sorted_states) > display_count:
        print(f"... ({len(sorted_states) - display_count} more patterns)")
    print()

    print("="*70)
    print("DENSITY OF STATES (ENERGY SPECTRUM)")
    print("="*70)
    print(f"{'Contacts':<12} {'Energy (kT)':<15} {'Conformations':>15} {'Fraction':>12}")
    print("-"*70)

    # Display density of states sorted by energy
    for c in sorted(contacts.keys()):
        energy = config.eps * c
        count = contacts[c]
        fraction = count / nconfs
        print(f"{c:<12} {energy:<15.1f} {count:>15} {fraction:>12.6f}")

    print("-"*70)
    print(f"{'TOTAL':<12} {'':<15} {nconfs:>15} {1.0:>12.6f}")
    print()
    print(f"Temperature: {config.T:.1f} K")
    print(f"Lattice: 2D square lattice")
    print("="*70)

    # ------------------------------------------------------------------------
    # 7. SUMMARY STATISTICS
    # ------------------------------------------------------------------------

    print("\n" + "="*70)
    print("SUMMARY STATISTICS - 2D SQUARE LATTICE")
    print("="*70)

    # Find ground state (most contacts = lowest energy)
    if contacts:
        max_contacts = max(contacts.keys())
        ground_state_energy = config.eps * max_contacts
        ground_state_degeneracy = contacts[max_contacts]

        print(f"Ground state energy: {ground_state_energy:.1f} kT ({max_contacts} contacts)")
        print(f"Ground state degeneracy: {ground_state_degeneracy}")

        import math
        print(f"Ground state entropy: kB * ln({ground_state_degeneracy}) = {math.log(ground_state_degeneracy):.4f} kB")

    # Find unfolded state (zero contacts = highest energy)
    if 0 in contacts:
        unfolded_degeneracy = contacts[0]
        print(f"Unfolded state (0 contacts) degeneracy: {unfolded_degeneracy}")
        print(f"Unfolded state entropy: kB * ln({unfolded_degeneracy}) = {math.log(unfolded_degeneracy):.4f} kB")

    # Total conformational entropy
    total_entropy = math.log(nconfs)
    print(f"Total conformational entropy: kB * ln({nconfs}) = {total_entropy:.4f} kB")

    print(f"\nNote: Results are for 2D square lattice with 4 possible move directions")
    print(f"Conformational space scaling: approximately 3^n = 3^{replicas[0].chain.n} configurations")

    print("="*70)
    print(f"\nResults written to: {combined_path}")
    print("Enumeration complete.\n")

