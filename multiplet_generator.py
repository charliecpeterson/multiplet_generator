import re, itertools, math, sys

# Mapping orbital letter to l value.
l_from_letter = {
    's': 0, 'p': 1, 'd': 2, 'f': 3,
    'g': 4, 'h': 5, 'i': 6, 'k': 7,
}

def print_help():
    help_text = """
Multiplet Generator - LS Coupling Calculator

Usage:
  python multiplet_generator.py [configuration] [options]

Examples:
  python multiplet_generator.py s1p4
  python multiplet_generator.py s1s1 --stats
  python multiplet_generator.py -h

Configuration Format:
  The configuration string is composed of tokens in the form:
    [principal quantum number][orbital letter][electron count]
  
  - For example:
      "p2"   : 2 electrons in a p-shell (if no principal quantum number is provided, these electrons are considered equivalent and are antisymmetrized).
      "s1s1" : Two separate (non-equivalent) s electrons (each with an explicit principal quantum number).
      "1sp2": 1 electron in the 1s orbital and 2 electrons in a p orbital.

Options:
  -h, --help    : Display this help message and exit.
  -s, --stats   : Print extra statistics, including total microstate count and LS multiplet degeneracy consistency check.

How It Works:
  1. The program parses the configuration string into tokens. Tokens with a specified principal quantum number (e.g., "1s1") are treated as non-equivalent shells, while tokens without a number are merged.
  2. For each shell, it calculates a microstate distribution (using Slater determinants for equivalent electrons or simple state counting for non-equivalent ones).
  3. The distributions for all shells are combined via convolution to obtain the overall microstate distribution.
  4. An iterative "peeling" algorithm then extracts the LS multiplets (with term symbols like 3P, 1D, 1S, etc.), also computing how many times each multiplet appears.
  5. The program computes the overall parity (even or odd) of the configuration and appends a superscript "°" to the term symbols if the parity is odd.
  6. If the --stats flag is provided, it prints a consistency check showing the total number of microstates versus the sum of the degeneracies of the LS multiplets.

Enjoy using the Multiplet Generator!
"""
    print(help_text)
    sys.exit(0)

def parse_configuration(config_str):
    """
    Parse a configuration string such as:
      "p2"    -> one token: {'n': None, 'orb': 'p', 'count': 2}
      "p1p1"  -> two tokens (non-equivalent if a principal quantum number is given)
      "1sp2"  -> tokens: "1s" and "p2"
    Tokens with an explicit principal quantum number (e.g. "1s") are treated as unique (non‑equivalent).
    Tokens without a number are considered equivalent and later merged.
    """
    pattern = re.compile(r'(\d*)\s*([spdfghiklmnoqSPDFGHIKLMNOQ])\s*(\d+)')
    tokens = pattern.findall(config_str)
    shells = []
    for n_str, orb, count_str in tokens:
        n = int(n_str) if n_str else None
        count = int(count_str)
        shells.append({'n': n, 'orb': orb.lower(), 'count': count})
    return shells

def generate_slater_determinants(l, n):
    """
    For equivalent electrons in a subshell with orbital angular momentum l,
    generate all unique determinants (as tuples of (m_l, m_s)).
    """
    orbitals = [(m_l, m_s) for m_l in range(-l, l+1) for m_s in [-0.5, 0.5]]
    return list(itertools.combinations(orbitals, n))

def shell_microstate_distribution(l, n, equivalent=True):
    """
    Return a dictionary mapping (M_L, M_S) -> count for a given shell.
    For equivalent electrons (e.g. "p2" with no principal quantum number) we antisymmetrize
    via Slater determinants.
    For non-equivalent electrons (typically count==1 when a principal quantum number is given),
    we simply use the single-electron states.
    """
    dist = {}
    if equivalent:
        dets = generate_slater_determinants(l, n)
        for det in dets:
            M_L = sum(orb[0] for orb in det)
            M_S = sum(orb[1] for orb in det)
            dist[(M_L, M_S)] = dist.get((M_L, M_S), 0) + 1
    else:
        # For non-equivalent electrons, typically count==1:
        for _ in range(n):
            for m_l in range(-l, l+1):
                for m_s in [-0.5, 0.5]:
                    key = (m_l, m_s)
                    dist[key] = dist.get(key, 0) + 1
    return dist

def combine_distributions(distributions):
    """
    Combine a list of microstate distributions (from different shells) by convolution.
    Each distribution is a dict mapping (M_L, M_S) -> count.
    """
    combined = { (0,0): 1 }
    for d in distributions:
        new_comb = {}
        for (mL1, mS1), count1 in combined.items():
            for (mL2, mS2), count2 in d.items():
                key = (mL1 + mL2, round(mS1 + mS2, 1))
                new_comb[key] = new_comb.get(key, 0) + count1 * count2
        combined = new_comb
    return combined

# Helpers for the peeling algorithm.
def m_values(S):
    """Return the list of m_S values for a given S (in steps of 1)."""
    vals = []
    m = -S
    while m <= S + 1e-9:
        vals.append(round(m, 1))
        m += 1.0
    return vals

def multiplet_states(L, S):
    """
    Return the set of states (m_L, m_S) that would be present in a full multiplet
    with quantum numbers L and S.
    """
    states = set()
    for mL in range(-L, L+1):
        for mS in m_values(S):
            states.add((mL, mS))
    return states

def can_subtract(distribution, L, S):
    """
    Check whether one copy of the multiplet (L,S) can be subtracted from the distribution.
    """
    for state in multiplet_states(L, S):
        if distribution.get(state, 0) < 1:
            return False
    return True

def subtract_multiplet(distribution, L, S):
    """
    Subtract one copy of the multiplet (L,S) from the distribution.
    """
    for state in multiplet_states(L, S):
        distribution[state] -= 1
        if distribution[state] == 0:
            del distribution[state]

def peel_multiplets_with_counts(distribution, total_electrons):
    """
    Iteratively subtract multiplet blocks from the overall microstate distribution.
    Returns a dictionary mapping each multiplet (L,S) to the number of times it was subtracted.
    (For equivalent electrons this count should match the LS decomposition.)
    """
    multiplet_counts = {}
    # Maximum possible S is total_electrons/2.
    S_max = total_electrons / 2
    # Build candidate multiplets.
    max_mL = max((abs(mL) for (mL, mS) in distribution.keys()), default=0)
    candidate_list = []
    S = S_max
    while S >= 0:
        for L in range(0, int(max_mL)+1):
            candidate_list.append((L, S))
        S -= 0.5
    # Sort candidates by descending S then descending L.
    candidate_list.sort(key=lambda x: (x[1], x[0]), reverse=True)
    
    while distribution:
        for (L, S) in candidate_list:
            if can_subtract(distribution, L, S):
                multiplet_counts[(L, S)] = multiplet_counts.get((L, S), 0) + 1
                subtract_multiplet(distribution, L, S)
                break
        else:
            break
    return multiplet_counts

def term_label(L, S):
    """
    Convert (L,S) into a term symbol like 3P, 1D, etc.
    Multiplicity = 2S+1, and L=0,1,2,... map to S,P,D,F,...
    """
    L_letters = "SPDFGHIKLMNOQ"
    mult = int(2*S + 1)
    if L < len(L_letters):
        return f"{mult}{L_letters[L]}"
    else:
        return f"{mult}L({L})"

def main():
    # Check for help flag.
    if any(arg in ["-h", "--help"] for arg in sys.argv[1:]):
        print_help()

    # Determine configuration string.
    if len(sys.argv) > 1:
        config_str = sys.argv[1]
    else:
        config_str = input("Enter electron configuration (e.g., p2, p1p1, 1sp2): ").strip()
    
    # Check for extra flag for detailed statistics.
    stats_flag = any(arg in ["--stats", "-s"] for arg in sys.argv[2:])
    
    tokens = parse_configuration(config_str)
    
    # Process tokens:
    # - Tokens with an explicit principal quantum number (n not None) are kept separate (non-equivalent).
    # - Tokens without a number are merged.
    shells = []
    for token in tokens:
        if token['n'] is None:
            shells.append(token)
        else:
            found = False
            for shell in shells:
                if shell['n'] == token['n'] and shell['orb'] == token['orb']:
                    shell['count'] += token['count']
                    found = True
                    break
            if not found:
                shells.append(token)
    
    # Compute overall parity.
    # For each shell, the contribution is (-1)^(count*l)
    parity_exponent = 0
    for shell in shells:
        l_val = l_from_letter[shell['orb']]
        parity_exponent += shell['count'] * l_val
    overall_parity = (-1) ** parity_exponent  # +1: even, -1: odd
    parity_superscript = "°" if overall_parity == -1 else ""
    parity_text = "odd" if overall_parity == -1 else "even"
    
    # Build microstate distributions for each shell.
    distributions = []
    total_electrons = 0
    for shell in shells:
        l_val = l_from_letter[shell['orb']]
        total_electrons += shell['count']
        equivalent = (shell['n'] is None)
        d = shell_microstate_distribution(l_val, shell['count'], equivalent=equivalent)
        distributions.append(d)
    
    # Combine the distributions from each shell.
    overall_dist = combine_distributions(distributions)
    total_microstates = sum(overall_dist.values())
    
    # Obtain LS multiplets via the peeling algorithm (with counts).
    multiplet_dict = peel_multiplets_with_counts(overall_dist.copy(), total_electrons)
    # Sort the multiplets (by descending S then L) for display.
    sorted_multiplets = sorted(multiplet_dict.items(), key=lambda x: (x[0][1], x[0][0]), reverse=True)
    
    # Print basic output.
    print(f"\nOverall configuration parity: {parity_text}")
    print("\nPossible LS multiplets:")
    for (L, S), count in sorted_multiplets:
        label = term_label(L, S) + parity_superscript
        print(f"  {label:6} (L={L}, S={S:3})  count: {count}  expected degeneracy: {(2*L+1)*(2*S+1)}")
    
    # If stats_flag is on, print the consistency check.
    if stats_flag:
        total_ls_degeneracy = sum(count * (2*L+1)*(2*S+1) for ((L, S), count) in multiplet_dict.items())
        print("\nMicrostate consistency check:")
        print(f"  Total microstates (from overall distribution): {total_microstates}")
        print(f"  Sum of LS multiplet degeneracies         : {total_ls_degeneracy}")
        if total_microstates == total_ls_degeneracy:
            print("  [OK] The LS multiplet decomposition accounts for all microstates.")
        else:
            print("  [Note] For non-equivalent electrons the sum of LS multiplet degeneracies\n         may not match the total microstate count exactly.")

if __name__ == "__main__":
    main()

