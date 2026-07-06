#!/usr/bin/env python3
# LS-coupling multiplet tool: terms, J levels, Hund's-rules ground state,
# microstate tables, and the classic term-energy expressions for p^n / d^n.
# A stand-in for Table 13-2 of Slater, "Quantum Theory of Atomic Structure".

import argparse
import itertools
import math
import re
import sys
from fractions import Fraction

L_OF = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6, 'k': 7}
# Letters for L = 0, 1, 2, ...; the series skips J and does not reuse S, P, D, F.
L_LETTERS = "SPDFGHIKLMNOQRTUVWXYZ"

TOKEN = re.compile(r'(\d*)([a-zA-Z])(\d+)')


def fail(msg):
    sys.exit(f"error: {msg}")


def capacity(l):
    return 4 * l + 2


def parse_config(chunks):
    """Parse tokens like p2, 3d4, 2p1 into subshells.

    Tokens with the same explicit n and letter merge into one subshell
    (2p1 2p1 == 2p2). Tokens without an n are each a distinct subshell,
    so 'p1 p1' means two nonequivalent p electrons.
    """
    shells = []

    def add(n, orb, count):
        l = L_OF[orb]
        if count > capacity(l):
            fail(f"'{orb}{count}': a {orb} subshell holds at most {capacity(l)} "
                 f"electrons. If you meant several subshells, separate the "
                 f"tokens (e.g. '2p1 3p1', not '2p13p1').")
        if n is not None and n <= l:
            fail(f"'{n}{orb}' does not exist (n must be at least l+1)")
        if n is not None:
            for s in shells:
                if s['n'] == n and s['orb'] == orb:
                    s['count'] += count
                    if s['count'] > capacity(l):
                        fail(f"{n}{orb} holds at most {capacity(l)} electrons "
                             f"({s['count']} given in total)")
                    return
        shells.append({'n': n, 'orb': orb, 'l': l, 'count': count})

    for chunk in re.split(r'[,./\s]+', ' '.join(chunks)):
        if not chunk:
            continue
        if not re.fullmatch(r'(?:\d*[a-zA-Z]\d+)+', chunk):
            fail(f"cannot parse '{chunk}' (tokens look like p2, 3d4, 2p1)")
        for n_str, orb, count_str in TOKEN.findall(chunk):
            orb = orb.lower()
            if orb not in L_OF:
                fail(f"unknown orbital letter '{orb}'")
            add(int(n_str) if n_str else None, orb, int(count_str))
    if not shells:
        fail("empty configuration")
    return shells


def shell_name(s):
    return f"{s['n'] if s['n'] is not None else ''}{s['orb']}{s['count']}"


def subshell_distribution(l, count):
    """Map (M_L, 2*M_S) -> number of Slater determinants for `count`
    equivalent electrons in a subshell of angular momentum l."""
    spin_orbitals = [(ml, ms2) for ml in range(-l, l + 1) for ms2 in (1, -1)]
    dist = {}
    for det in itertools.combinations(spin_orbitals, count):
        key = (sum(ml for ml, _ in det), sum(ms2 for _, ms2 in det))
        dist[key] = dist.get(key, 0) + 1
    return dist


def convolve(dists):
    total = {(0, 0): 1}
    for d in dists:
        merged = {}
        for (ml1, ms1), c1 in total.items():
            for (ml2, ms2), c2 in d.items():
                key = (ml1 + ml2, ms1 + ms2)
                merged[key] = merged.get(key, 0) + c1 * c2
        total = merged
    return total


def extract_terms(dist):
    """Decompose an (M_L, 2*M_S) distribution into LS terms.

    The number of terms with quantum numbers (L, S) is
    n(L,S) - n(L+1,S) - n(L,S+1) + n(L+1,S+1), where n(ML,MS) is the
    microstate count: Slater's Sec. 13-2 bookkeeping in closed form.
    """
    terms = {}
    max_l = max(ml for ml, _ in dist)
    max_s2 = max(ms2 for _, ms2 in dist)
    for L in range(max_l + 1):
        for s2 in range(max_s2 % 2, max_s2 + 1, 2):
            c = (dist.get((L, s2), 0) - dist.get((L + 1, s2), 0)
                 - dist.get((L, s2 + 2), 0) + dist.get((L + 1, s2 + 2), 0))
            if c < 0:
                raise AssertionError("inconsistent microstate distribution")
            if c:
                terms[(L, s2)] = c
    return terms


def half(two_x):
    return str(Fraction(two_x, 2))


def plural(count, noun):
    return f"{count} {noun}{'' if count == 1 else 's'}"


def term_label(L, s2, odd=False):
    letter = L_LETTERS[L] if L < len(L_LETTERS) else f"(L={L})"
    return f"{s2 + 1}{letter}" + ("°" if odd else "")


def j_values(L, s2):
    """Allowed 2*J for a term, from |L-S| to L+S."""
    return list(range(abs(2 * L - s2), 2 * L + s2 + 1, 2))


def lande_g(L, s2, j2):
    if j2 == 0:
        return None
    J, S = j2 / 2, s2 / 2
    return 1 + (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1))


def open_subshells(shells):
    return [s for s in shells if 0 < s['count'] < capacity(s['l'])]


def print_terms(terms, odd, total_microstates):
    print("\nTerms:")
    ordered = sorted(terms.items(), key=lambda kv: (-kv[0][1], -kv[0][0]))
    width = max(len(term_label(L, s2, odd)) for (L, s2) in terms)
    for (L, s2), count in ordered:
        label = term_label(L, s2, odd)
        js = ", ".join(half(j2) for j2 in j_values(L, s2))
        states = (2 * L + 1) * (s2 + 1) * count
        mult = f"×{count}" if count > 1 else "  "
        print(f"  {label:<{width}} {mult}   J = {js:<24} "
              f"({plural(states, 'state')})")
    n_terms = sum(terms.values())
    print(f"  -- {plural(n_terms, 'term')}, "
          f"{plural(total_microstates, 'microstate')}")


def print_hund(shells, terms, odd):
    opens = open_subshells(shells)
    if not opens:
        print("\nGround state: 1S0 (all subshells closed)")
        return
    if len(opens) > 1:
        print("\nHund's rules apply to a single open subshell; "
              "skipping ground-state analysis.")
        return
    l, n = opens[0]['l'], opens[0]['count']
    L, s2 = max(terms, key=lambda t: (t[1], t[0]))
    half_fill = 2 * l + 1
    j2 = abs(2 * L - s2) if n <= half_fill else 2 * L + s2
    label = term_label(L, s2, odd)
    if n < half_fill:
        note = "normal multiplet (lowest J lies lowest): shell less than half-filled"
    elif n > half_fill:
        note = "inverted multiplet (highest J lies lowest): shell more than half-filled"
    else:
        note = "half-filled shell: L = 0, single J value"
    print(f"\nGround state (Hund's rules): {label}, J = {half(j2)}")
    print(f"  {note}")
    gs = [(j2v, lande_g(L, s2, j2v)) for j2v in j_values(L, s2)]
    parts = [f"g({half(j2v)}) = {g:.3f}" for j2v, g in gs if g is not None]
    if parts:
        print(f"  Landé g of the {label} levels: " + ", ".join(parts))
    S = s2 / 2
    mu = f"  μ_eff: spin-only 2√(S(S+1)) = {2 * math.sqrt(S * (S + 1)):.2f} μ_B"
    g0 = lande_g(L, s2, j2)
    if g0 is not None and L > 0:
        J = j2 / 2
        mu += f";  free-ion g_J√(J(J+1)) = {g0 * math.sqrt(J * (J + 1)):.2f} μ_B"
    elif j2 == 0:
        mu += ";  free-ion J = 0, μ_eff = 0"
    print(mu)


def couple(L1, s21, L2, s22):
    """Clebsch-Gordan series on L and S (spins doubled)."""
    return [(L, s2) for L in range(abs(L1 - L2), L1 + L2 + 1)
            for s2 in range(abs(s21 - s22), s21 + s22 + 1, 2)]


def sorted_terms(terms):
    return sorted(terms.items(), key=lambda kv: (-kv[0][1], -kv[0][0]))


def print_parentage(shells, terms, odd):
    """Couple the open subshells' own terms left to right, so every final
    term carries its parentage — the (3F)4F-style labels GRASP/NIST use."""
    opens = open_subshells(shells)
    if len(opens) < 2:
        return

    def sub_terms(s):
        return sorted_terms(
            extract_terms(subshell_distribution(s['l'], s['count'])))

    tally = {}

    def walk(parent, pcount, plabel, rest, indent):
        s, last = rest[0], len(rest) == 1
        for (L2, s22), c2 in sub_terms(s):
            kids = sorted(couple(parent[0], parent[1], L2, s22),
                          key=lambda t: (-t[1], -t[0]))
            shown = " ".join(term_label(L, s2, odd and last) for L, s2 in kids)
            mult = f"   ×{pcount * c2}" if pcount * c2 > 1 else ""
            print(f"{' ' * indent}{plabel} × "
                  f"{shell_name(s)}({term_label(L2, s22)}) → {shown}{mult}")
            for kid in kids:
                if last:
                    tally[kid] = tally.get(kid, 0) + pcount * c2
                else:
                    walk(kid, pcount * c2, f"({term_label(*kid)})",
                         rest[1:], indent + 2)

    print("\nParentage (subshell terms coupled left to right):")
    first = opens[0]
    for (L, s2), c in sub_terms(first):
        walk((L, s2), c, f"{shell_name(first)}({term_label(L, s2)})",
             opens[1:], 2)
    assert tally == dict(terms), "parentage coupling disagrees with terms"


def print_interval(shells, terms, odd, zeta):
    """Fine structure of the ground term from the Landé interval rule,
    E(J) - E(J-1) = λJ with λ = ±ζ/2S for a single open subshell."""
    opens = open_subshells(shells)
    if len(opens) != 1:
        print("\nInterval rule needs a single open subshell.")
        return
    l, n = opens[0]['l'], opens[0]['count']
    L, s2 = max(terms, key=lambda t: (t[1], t[0]))
    label = term_label(L, s2, odd)
    if L == 0 or s2 == 0:
        print(f"\n{label}: no first-order spin-orbit splitting "
              "(L = 0 or S = 0).")
        return
    lam = zeta / s2 * (1 if n < 2 * l + 1 else -1)
    levels = sorted((lam / 2 * j2 * (j2 + 2) / 4, j2) for j2 in j_values(L, s2))
    e0 = levels[0][0]
    print(f"\nInterval-rule fine structure of {label} "
          f"(ζ = {zeta:g}, λ = ±ζ/2S = {lam:+g}):")
    prev = None
    for e, j2 in levels:
        gap = f"   (gap {e - prev:.1f})" if prev is not None else ""
        print(f"  {label}{half(j2)}   E = {e - e0:10.1f}{gap}")
        prev = e


def jj_shell_options(l, n):
    """All (n_minus, n_plus) splits of n electrons over j = l-1/2, l+1/2."""
    opts = []
    for n_m in range(min(n, 2 * l), -1, -1):
        n_p = n - n_m
        if 0 <= n_p <= 2 * l + 2:
            opts.append((n_m, n_p))
    return opts


def jshell_m_dist(j2, count):
    """{2*M_J: count} for `count` equivalent electrons in a j-shell."""
    dist = {}
    for det in itertools.combinations(range(-j2, j2 + 1, 2), count):
        m = sum(det)
        dist[m] = dist.get(m, 0) + 1
    return dist


def extract_j_levels(dists):
    """Convolve {2M: count} distributions and peel into {2J: levels}."""
    total = {0: 1}
    for d in dists:
        merged = {}
        for m1, c1 in total.items():
            for m2, c2 in d.items():
                merged[m1 + m2] = merged.get(m1 + m2, 0) + c1 * c2
        total = merged
    levels = {}
    top = max(total)
    for j2 in range(top % 2, top + 1, 2):
        c = total.get(j2, 0) - total.get(j2 + 2, 0)
        if c < 0:
            raise AssertionError("inconsistent M_J distribution")
        if c:
            levels[j2] = c
    return levels


def print_jj(shells, terms, odd):
    """Relativistic (jj) subshell splits with their J values, plus a
    levels-per-J cross-check against the LS terms. The per-J totals are
    what GRASP's rcsfgenerate/rlevels report for the same configuration."""
    print("\njj coupling (x- means j = l-1/2, x means j = l+1/2):")
    opens = open_subshells(shells)
    if not opens:
        print("  closed shells only: a single CSF with J = 0")
        return
    per_shell = []
    for s in opens:
        name = f"{s['n'] if s['n'] is not None else ''}{s['orb']}"
        opts = []
        for n_m, n_p in jj_shell_options(s['l'], s['count']):
            parts, occ = [], []
            if n_m:
                parts.append(f"({name}-)^{n_m}")
                occ.append((2 * s['l'] - 1, n_m))
            if n_p:
                parts.append(f"({name})^{n_p}")
                occ.append((2 * s['l'] + 1, n_p))
            opts.append((" ".join(parts), occ))
        per_shell.append(opts)
    rows = []
    census = {}
    for combo in itertools.product(*per_shell):
        label = "  ".join(lbl for lbl, _ in combo)
        dists = [jshell_m_dist(j2, c) for _, occ in combo for j2, c in occ]
        levels = extract_j_levels(dists)
        for j2, c in levels.items():
            census[j2] = census.get(j2, 0) + c
        js = ", ".join(half(j2) + (f"×{c}" if c > 1 else "")
                       for j2, c in sorted(levels.items()))
        rows.append((label, js))
    width = max(len(label) for label, _ in rows)
    for label, js in rows:
        print(f"  {label:<{width}}   J = {js}")
    ls_census = {}
    for (L, s2), c in terms.items():
        for j2 in j_values(L, s2):
            ls_census[j2] = ls_census.get(j2, 0) + c
    fmt = lambda cen: "  ".join(f"J={half(j2)}: {c}"
                                for j2, c in sorted(cen.items()))
    print(f"  Levels per J (jj): {fmt(census)}   "
          f"({plural(sum(census.values()), 'level')})")
    verdict = "consistent" if census == ls_census else "MISMATCH — bug!"
    print(f"  Levels per J (LS): {fmt(ls_census)}   {verdict}")


# Reduction of the rotation-group character chi_L onto the irreps of O;
# with the g/u label from the configuration parity this is the weak-field
# splitting of each free-ion term in Oh.
O_CLASSES = ((1, 0.0), (8, 2 * math.pi / 3), (3, math.pi),
             (6, math.pi / 2), (6, math.pi))
O_IRREPS = (("A1", (1, 1, 1, 1, 1)), ("A2", (1, 1, 1, -1, -1)),
            ("E", (2, -1, 2, 0, 0)), ("T1", (3, 0, -1, 1, -1)),
            ("T2", (3, 0, -1, -1, 1)))


def oh_split(L):
    def chi(theta):
        if theta == 0.0:
            return float(2 * L + 1)
        return math.sin((2 * L + 1) * theta / 2) / math.sin(theta / 2)

    chars = [chi(theta) for _, theta in O_CLASSES]
    out = []
    for name, row in O_IRREPS:
        raw = sum(g * r * c
                  for (g, _), r, c in zip(O_CLASSES, row, chars)) / 24
        n = round(raw)
        if abs(raw - n) > 1e-6:
            raise AssertionError(f"non-integer reduction for L={L}")
        if n:
            out.append((name, n))
    dims = {'A1': 1, 'A2': 1, 'E': 2, 'T1': 3, 'T2': 3}
    assert sum(dims[name] * n for name, n in out) == 2 * L + 1
    return out


def print_oh(terms, odd):
    suffix = "u" if odd else "g"
    print("\nWeak-field splitting in Oh:")
    for (L, s2), count in sorted_terms(terms):
        pieces = " + ".join(
            f"{s2 + 1}{name}{suffix}" + (f"(×{n})" if n > 1 else "")
            for name, n in oh_split(L))
        mult = f"   ×{count}" if count > 1 else ""
        print(f"  {term_label(L, s2, odd):<5} → {pieces}{mult}")


def print_grid(dist):
    """M_L rows vs M_S columns, like Slater's Figs. 13-3/13-4."""
    mls = sorted({ml for ml, _ in dist}, reverse=True)
    ms2s = sorted({ms2 for _, ms2 in dist})
    headers = [half(ms2) for ms2 in ms2s]
    col = max(3, max(len(h) for h in headers) + 1)
    row_label = max(len(str(ml)) for ml in mls) + 5
    print("\nMicrostates by (M_L, M_S):")
    print(" " * row_label + "M_S " + "".join(f"{h:>{col}}" for h in headers))
    for ml in mls:
        cells = "".join(
            f"{dist.get((ml, ms2), '') or '':>{col}}" for ms2 in ms2s)
        print(f"{'M_L=' + str(ml):>{row_label}}    " + cells)


# Term energies from the diagonal sum rule (Slater Secs. 13-1..13-3; the
# d^n forms are the standard Racah-parameter tables, e.g. Griffith,
# "The Theory of Transition-Metal Ions", Table 4.6). A contribution common
# to every term of the configuration (F0 terms / multiples of Racah A) is
# omitted: it shifts all terms equally and cancels in splittings.
ENERGY_TABLES = {
    (1, 2): {
        'params': ('F2',),
        'terms': [
            ("3P", "-5 F2", lambda p: -5 * p['F2']),
            ("1D", "F2", lambda p: p['F2']),
            ("1S", "10 F2", lambda p: 10 * p['F2']),
        ],
        'note': "F2 = F^2/25 (Condon-Shortley); common F0 omitted.",
    },
    (1, 3): {
        'params': ('F2',),
        'terms': [
            ("4S", "-15 F2", lambda p: -15 * p['F2']),
            ("2D", "-6 F2", lambda p: -6 * p['F2']),
            ("2P", "0", lambda p: 0.0),
        ],
        'note': "F2 = F^2/25 (Condon-Shortley); common 3 F0 omitted.",
    },
    (2, 2): {
        'params': ('B', 'C'),
        'terms': [
            ("3F", "-8B", lambda p: -8 * p['B']),
            ("1D", "-3B + 2C", lambda p: -3 * p['B'] + 2 * p['C']),
            ("3P", "7B", lambda p: 7 * p['B']),
            ("1G", "4B + 2C", lambda p: 4 * p['B'] + 2 * p['C']),
            ("1S", "14B + 7C", lambda p: 14 * p['B'] + 7 * p['C']),
        ],
        'note': "Racah parameters; common A omitted.",
    },
    (2, 3): {
        'params': ('B', 'C'),
        'terms': [
            ("4F", "-15B", lambda p: -15 * p['B']),
            ("4P", "0", lambda p: 0.0),
            ("2G", "-11B + 3C", lambda p: -11 * p['B'] + 3 * p['C']),
            ("2P", "-6B + 3C", lambda p: -6 * p['B'] + 3 * p['C']),
            ("2H", "-6B + 3C", lambda p: -6 * p['B'] + 3 * p['C']),
            ("2F", "9B + 3C", lambda p: 9 * p['B'] + 3 * p['C']),
            ("2D", "5B + 5C - sqrt(193B^2 + 8BC + 4C^2)",
             lambda p: 5 * p['B'] + 5 * p['C'] - _d3_root(p)),
            ("2D", "5B + 5C + sqrt(193B^2 + 8BC + 4C^2)",
             lambda p: 5 * p['B'] + 5 * p['C'] + _d3_root(p)),
        ],
        'note': ("Racah parameters; common 3A omitted. "
                 "2P and 2H are degenerate at this level of treatment."),
    },
}


def _d3_root(p):
    B, C = p['B'], p['C']
    return math.sqrt(193 * B * B + 8 * B * C + 4 * C * C)


def print_energies(shells, odd, params):
    opens = open_subshells(shells)
    if len(opens) != 1:
        print("\nTerm energies are tabulated only for a single open subshell.")
        return
    l, n = opens[0]['l'], opens[0]['count']
    n_eff = min(n, capacity(l) - n)
    orb = opens[0]['orb']
    if n_eff <= 1:
        print(f"\n{orb}{n} is a single electron or hole: no electron-repulsion "
              "splitting, the configuration is one term.")
        return
    table = ENERGY_TABLES.get((l, n_eff))
    if table is None:
        print(f"\nNo energy table built in for {orb}{n_eff} "
              "(only p2, p3, d2, d3 and their hole partners).")
        return
    title = f"\nTerm energies for {orb}{n}"
    if n_eff != n:
        title += f" (hole equivalence: same splittings as {orb}{n_eff})"
    print(title + ":")
    have = all(params.get(name) is not None for name in table['params'])
    rows = []
    for label, expr, fn in table['terms']:
        value = fn(params) if have else None
        rows.append((label + ("°" if odd else ""), expr, value))
    if have:
        e0 = min(v for _, _, v in rows)
        rows.sort(key=lambda r: r[2])
    width = max(len(r[0]) for r in rows)
    expr_w = max(len(r[1]) for r in rows)
    for label, expr, value in rows:
        line = f"  {label:<{width}}   {expr:<{expr_w}}"
        if have:
            line += f"   E - E0 = {value - e0:12.1f}"
        print(line)
    print(f"  ({table['note']})")
    if not have:
        needed = " and ".join(f"--{name}" for name in table['params'])
        print(f"  (pass {needed} to evaluate numerically, "
              "in your unit of choice, e.g. cm^-1)")


def main():
    ap = argparse.ArgumentParser(
        description="LS terms, J levels, and Hund's-rules ground state for an "
                    "electron configuration (Slater Table 13-2 and friends).",
        epilog="Configuration tokens: [n]<letter><count>, e.g. 'p2', '3d4', "
               "'2p1 3p1'. Same n+letter merge into one subshell; tokens "
               "without n are distinct subshells. Electrons within a token "
               "are equivalent (Pauli applies); different tokens couple as "
               "independent shells.")
    ap.add_argument('config', nargs='+',
                    help="configuration, e.g. p2 | 3d4 | '2p1 3p1' | 1s2.2p3")
    ap.add_argument('--grid', action='store_true',
                    help='print the M_L/M_S microstate table')
    ap.add_argument('--energies', action='store_true',
                    help='print term-energy expressions (p2/p3 in F2, d2/d3 '
                         'in Racah B,C; hole partners map onto these)')
    ap.add_argument('--B', type=float, help='Racah B (implies --energies)')
    ap.add_argument('--C', type=float, help='Racah C (implies --energies)')
    ap.add_argument('--F2', type=float,
                    help='Slater-Condon F2 for p^n (implies --energies)')
    ap.add_argument('--jj', action='store_true',
                    help='relativistic (jj) subshell splits with J values '
                         'and a per-J cross-check against LS — matches '
                         "GRASP's CSF/level bookkeeping")
    ap.add_argument('--zeta', type=float, metavar='Z',
                    help='spin-orbit ζ: interval-rule fine structure of the '
                         'ground term (same unit as ζ)')
    ap.add_argument('--oh', action='store_true',
                    help='weak-field splitting of each term in Oh')
    args = ap.parse_args()

    shells = parse_config(args.config)
    odd = sum(s['l'] * s['count'] for s in shells) % 2 == 1
    dist = convolve([subshell_distribution(s['l'], s['count']) for s in shells])
    terms = extract_terms(dist)
    total = sum(dist.values())
    assert total == sum(c * (2 * L + 1) * (s2 + 1)
                        for (L, s2), c in terms.items())

    names = " ".join(shell_name(s) for s in shells)
    print(f"Configuration: {names}   "
          f"({plural(total, 'microstate')}, parity {'odd' if odd else 'even'})")
    for s in open_subshells(shells):
        n_holes = capacity(s['l']) - s['count']
        if s['count'] > capacity(s['l']) // 2 and n_holes > 0:
            print(f"  {s['orb']}{s['count']} has the same terms as "
                  f"{s['orb']}{n_holes} (hole equivalence)")

    print_terms(terms, odd, total)
    print_hund(shells, terms, odd)
    print_parentage(shells, terms, odd)
    if args.grid:
        print_grid(dist)
    if args.zeta is not None:
        print_interval(shells, terms, odd, args.zeta)
    if args.jj:
        print_jj(shells, terms, odd)
    if args.oh:
        print_oh(terms, odd)
    if args.energies or any(v is not None for v in (args.B, args.C, args.F2)):
        print_energies(shells, odd,
                       {'B': args.B, 'C': args.C, 'F2': args.F2})


if __name__ == "__main__":
    main()
