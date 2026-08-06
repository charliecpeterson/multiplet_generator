"""Enumerate LS terms and J-parity levels for atomic configurations.

The CLI is intended as a symmetry and state-count sanity check before running
an atomic-structure code. It deliberately does not predict level energies.
"""

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import urlencode

from grasp_preflight import (
    JjAnalysis,
    ParentagePath,
    TransitionResult,
    allowed_e1_transitions,
    j_shell_m_distribution,
    jj_analysis,
    parentage_paths,
)


ORBITAL_LETTERS = "spdfghiklmnoq"
L_FROM_LETTER = {letter: l_value for l_value, letter in enumerate(ORBITAL_LETTERS)}
TERM_LETTERS = "SPDFGHIKLMNOQRTUVWXYZ"
TOKEN_PATTERN = re.compile(
    rf"(?P<n>\d*)(?P<orb>[{ORBITAL_LETTERS}])(?:\^)?(?P<count>\d+)",
    re.IGNORECASE,
)

Weight = Tuple[int, int]
Distribution = Dict[Weight, int]
SubshellKey = Tuple[Optional[int], int]


class ConfigurationError(ValueError):
    """Raised when a configuration cannot represent valid atomic subshells."""


@dataclass(frozen=True)
class Subshell:
    n: Optional[int]
    l: int
    occupancy: int

    @property
    def capacity(self) -> int:
        return 2 * (2 * self.l + 1)

    @property
    def effective_occupancy(self) -> int:
        return min(self.occupancy, self.capacity - self.occupancy)

    @property
    def is_hole_like(self) -> bool:
        return self.occupancy > self.capacity // 2

    @property
    def is_open(self) -> bool:
        return 0 < self.occupancy < self.capacity

    @property
    def name(self) -> str:
        principal = "" if self.n is None else str(self.n)
        return f"{principal}{ORBITAL_LETTERS[self.l]}"

    @property
    def label(self) -> str:
        return f"{self.name}{self.occupancy}"


@dataclass(frozen=True)
class Configuration:
    subshells: Tuple[Subshell, ...]

    @property
    def electron_count(self) -> int:
        return sum(shell.occupancy for shell in self.subshells)

    @property
    def parity(self) -> int:
        exponent = sum(shell.l * shell.occupancy for shell in self.subshells)
        return -1 if exponent % 2 else 1

    @property
    def parity_name(self) -> str:
        return "odd" if self.parity == -1 else "even"

    @property
    def parity_symbol(self) -> str:
        return "°" if self.parity == -1 else ""

    @property
    def label(self) -> str:
        return " ".join(shell.label for shell in self.subshells)


@dataclass(frozen=True)
class Term:
    L: int
    two_s: int
    occurrences: int

    @property
    def multiplicity(self) -> int:
        return self.two_s + 1

    @property
    def label(self) -> str:
        if self.L < len(TERM_LETTERS):
            orbital = TERM_LETTERS[self.L]
        else:
            orbital = f"L({self.L})"
        return f"{self.multiplicity}{orbital}"

    @property
    def degeneracy(self) -> int:
        return (2 * self.L + 1) * (self.two_s + 1)


@dataclass(frozen=True)
class HundResult:
    term: Term
    two_j: int
    filling: str


def parse_configuration(config_text: str) -> Configuration:
    """Parse compact forms such as ``2p2 3s1`` or the shorthand ``p2s1``."""
    if not config_text.strip():
        raise ConfigurationError("configuration is empty")

    parsed: List[Tuple[Optional[int], int, int]] = []
    cursor = 0
    for match in TOKEN_PATTERN.finditer(config_text):
        gap = config_text[cursor : match.start()]
        if gap.strip(" ,"):
            raise ConfigurationError(f"cannot parse configuration near {gap!r}")

        n_text = match.group("n")
        n = int(n_text) if n_text else None
        orbital = match.group("orb").lower()
        l_value = L_FROM_LETTER[orbital]
        occupancy = int(match.group("count"))

        if occupancy == 0:
            raise ConfigurationError("subshell occupancies must be positive")
        if n is not None and l_value >= n:
            raise ConfigurationError(
                f"{n}{orbital} is invalid: an n={n} shell requires l < {n}"
            )

        parsed.append((n, l_value, occupancy))
        cursor = match.end()

    trailing = config_text[cursor:]
    if trailing.strip(" ,"):
        raise ConfigurationError(f"cannot parse configuration near {trailing!r}")
    if not parsed:
        raise ConfigurationError(f"cannot parse configuration {config_text!r}")

    occupancies: Dict[SubshellKey, int] = {}
    order: List[SubshellKey] = []
    for n, l_value, occupancy in parsed:
        key = (n, l_value)
        if key not in occupancies:
            occupancies[key] = 0
            order.append(key)
        occupancies[key] += occupancy

    subshells = []
    for n, l_value in order:
        occupancy = occupancies[(n, l_value)]
        capacity = 2 * (2 * l_value + 1)
        if occupancy > capacity:
            label = f"{'' if n is None else n}{ORBITAL_LETTERS[l_value]}"
            raise ConfigurationError(
                f"{label}{occupancy} exceeds the subshell capacity of {capacity}"
            )
        subshells.append(Subshell(n=n, l=l_value, occupancy=occupancy))

    return Configuration(tuple(subshells))


def shell_microstate_distribution(l_value: int, occupancy: int) -> Distribution:
    """Count determinant weights without materializing the determinants themselves."""
    capacity = 2 * (2 * l_value + 1)
    if not 0 <= occupancy <= capacity:
        raise ValueError(f"occupancy must be between 0 and {capacity}")

    effective_occupancy = min(occupancy, capacity - occupancy)
    states: Dict[Tuple[int, int, int], int] = {(0, 0, 0): 1}
    for m_l in range(-l_value, l_value + 1):
        for two_m_s in (-1, 1):
            updated = dict(states)
            for (chosen, total_m_l, total_two_m_s), count in states.items():
                if chosen == effective_occupancy:
                    continue
                key = (chosen + 1, total_m_l + m_l, total_two_m_s + two_m_s)
                updated[key] = updated.get(key, 0) + count
            states = updated

    distribution: Distribution = {}
    for (chosen, total_m_l, total_two_m_s), count in states.items():
        if chosen == effective_occupancy:
            distribution[(total_m_l, total_two_m_s)] = count
    return distribution


def combine_distributions(distributions: Iterable[Mapping[Weight, int]]) -> Distribution:
    combined: Distribution = {(0, 0): 1}
    for distribution in distributions:
        updated: Distribution = {}
        for (m_l_left, two_m_s_left), left_count in combined.items():
            for (m_l_right, two_m_s_right), right_count in distribution.items():
                key = (m_l_left + m_l_right, two_m_s_left + two_m_s_right)
                updated[key] = updated.get(key, 0) + left_count * right_count
        combined = updated
    return combined


def configuration_microstates(configuration: Configuration) -> Distribution:
    return combine_distributions(
        shell_microstate_distribution(shell.l, shell.occupancy)
        for shell in configuration.subshells
    )


def decompose_ls(distribution: Mapping[Weight, int]) -> Tuple[Term, ...]:
    """Decompose an ML/MS weight distribution using highest-weight multiplicities."""
    if not distribution:
        raise ValueError("microstate distribution is empty")

    max_l = max(abs(m_l) for m_l, _ in distribution)
    max_two_s = max(abs(two_m_s) for _, two_m_s in distribution)
    spin_parity = next(iter(distribution))[1] % 2
    terms = []

    def weight(m_l: int, two_m_s: int) -> int:
        return distribution.get((m_l, two_m_s), 0)

    for two_s in range(spin_parity, max_two_s + 1, 2):
        for L in range(max_l + 1):
            occurrences = (
                weight(L, two_s)
                - weight(L + 1, two_s)
                - weight(L, two_s + 2)
                + weight(L + 1, two_s + 2)
            )
            if occurrences < 0:
                raise RuntimeError(
                    f"invalid LS decomposition at L={L}, 2S={two_s}: {occurrences}"
                )
            if occurrences:
                terms.append(Term(L=L, two_s=two_s, occurrences=occurrences))

    terms.sort(key=lambda term: (term.two_s, term.L), reverse=True)
    microstates = sum(distribution.values())
    term_states = sum(term.occurrences * term.degeneracy for term in terms)
    if term_states != microstates:
        raise RuntimeError(
            f"LS decomposition accounts for {term_states} of {microstates} microstates"
        )
    return tuple(terms)


def terms_for_configuration(configuration: Configuration) -> Tuple[Term, ...]:
    return decompose_ls(configuration_microstates(configuration))


def j_values(term: Term) -> Tuple[int, ...]:
    lower = abs(2 * term.L - term.two_s)
    upper = 2 * term.L + term.two_s
    return tuple(range(lower, upper + 1, 2))


def lande_g(term: Term, two_j: int) -> Optional[float]:
    if two_j == 0:
        return None
    J = two_j / 2
    S = term.two_s / 2
    return 1 + (
        J * (J + 1) + S * (S + 1) - term.L * (term.L + 1)
    ) / (2 * J * (J + 1))


def j_block_counts(terms: Sequence[Term]) -> Dict[int, int]:
    counts: Dict[int, int] = {}
    for term in terms:
        for two_j in j_values(term):
            counts[two_j] = counts.get(two_j, 0) + term.occurrences
    return dict(sorted(counts.items()))


def hund_ground(
    configuration: Configuration, terms: Sequence[Term]
) -> Tuple[Optional[HundResult], Optional[str]]:
    open_subshells = [shell for shell in configuration.subshells if shell.is_open]
    if not open_subshells:
        singlet_s = next(
            (term for term in terms if term.L == 0 and term.two_s == 0), None
        )
        if singlet_s is None:
            return None, "closed-shell configuration did not produce a 1S term"
        return HundResult(singlet_s, 0, "closed"), None
    if len(open_subshells) != 1:
        return None, "Hund ground-level guidance is limited to one open subshell"

    shell = open_subshells[0]
    maximum_spin = max(term.two_s for term in terms)
    maximum_l = max(term.L for term in terms if term.two_s == maximum_spin)
    ground_term = next(
        term
        for term in terms
        if term.two_s == maximum_spin and term.L == maximum_l
    )

    half = shell.capacity // 2
    if shell.occupancy < half:
        filling = "less than half-filled"
        two_j = abs(2 * ground_term.L - ground_term.two_s)
    elif shell.occupancy > half:
        filling = "more than half-filled"
        two_j = 2 * ground_term.L + ground_term.two_s
    else:
        filling = "half-filled"
        two_j = ground_term.two_s
    return HundResult(ground_term, two_j, filling), None


def format_half(two_value: int) -> str:
    if two_value % 2 == 0:
        return str(two_value // 2)
    return f"{two_value}/2"


def nist_levels_url(spectrum: str) -> str:
    query = urlencode(
        {
            "biblio": 1,
            "conf_out": 1,
            "j_out": 1,
            "lande_out": 1,
            "level_out": 1,
            "perc_out": 1,
            "spectrum": spectrum,
            "term_out": 1,
        }
    )
    return f"https://physics.nist.gov/cgi-bin/ASD/energy1.pl?{query}"


def analysis_dict(configuration: Configuration, terms: Sequence[Term]) -> dict:
    distribution = configuration_microstates(configuration)
    hund, hund_note = hund_ground(configuration, terms)
    term_entries = []
    for term in terms:
        levels = []
        for two_j in j_values(term):
            levels.append(
                {
                    "J": two_j / 2,
                    "degeneracy": two_j + 1,
                    "lande_g": lande_g(term, two_j),
                }
            )
        term_entries.append(
            {
                "term": term.label + configuration.parity_symbol,
                "L": term.L,
                "S": term.two_s / 2,
                "multiplicity": term.multiplicity,
                "occurrences": term.occurrences,
                "degeneracy_per_occurrence": term.degeneracy,
                "levels": levels,
            }
        )

    hund_entry = None
    if hund is not None:
        hund_entry = {
            "term": hund.term.label + configuration.parity_symbol,
            "J": hund.two_j / 2,
            "filling": hund.filling,
        }

    blocks = [
        {
            "J": two_j / 2,
            "parity": configuration.parity,
            "levels": count,
            "magnetic_sublevels": count * (two_j + 1),
        }
        for two_j, count in j_block_counts(terms).items()
    ]
    return {
        "configuration": configuration.label,
        "electron_count": configuration.electron_count,
        "parity": configuration.parity_name,
        "microstates": sum(distribution.values()),
        "terms": term_entries,
        "j_parity_blocks": blocks,
        "hund_ground": hund_entry,
        "hund_note": hund_note,
    }


def transition_dict(result: TransitionResult, initial_parity: int, final_parity: int) -> dict:
    entries = []
    initial_symbol = "°" if initial_parity == -1 else ""
    final_symbol = "°" if final_parity == -1 else ""
    for transition in result.transitions:
        entries.append(
            {
                "initial": (
                    f"{transition.initial_term.label}{initial_symbol}_"
                    f"{format_half(transition.initial_two_j)}"
                ),
                "final": (
                    f"{transition.final_term.label}{final_symbol}_"
                    f"{format_half(transition.final_two_j)}"
                ),
                "term_occurrence_combinations": transition.combinations,
            }
        )
    return {"allowed": entries, "reason": result.reason}


def terms_for_subshell(shell: Subshell) -> Tuple[Term, ...]:
    return terms_for_configuration(Configuration((shell,)))


def parentage_dict(paths: Sequence[ParentagePath], parity: int) -> list:
    parity_symbol = "°" if parity == -1 else ""
    return [
        {
            "parents": list(path.parents),
            "coupling_sequence": list(path.intermediate_terms),
            "final_term": path.final_label + parity_symbol,
            "occurrences": path.occurrences,
        }
        for path in paths
    ]


def jj_dict(analysis: JjAnalysis, parity: int) -> dict:
    return {
        "parity": "odd" if parity == -1 else "even",
        "relativistic_configurations": [
            {
                "configuration": row.label,
                "microstates": row.microstates,
                "J_levels": [
                    {"J": two_j / 2, "csfs": count}
                    for two_j, count in row.levels
                ],
            }
            for row in analysis.rows
        ],
        "jj_census": [
            {"J": two_j / 2, "levels": count}
            for two_j, count in analysis.jj_census
        ],
        "ls_census": [
            {"J": two_j / 2, "levels": count}
            for two_j, count in analysis.ls_census
        ],
        "consistent": analysis.consistent,
    }


def print_analysis(configuration: Configuration, terms: Sequence[Term], stats: bool) -> None:
    distribution = configuration_microstates(configuration)
    microstates = sum(distribution.values())
    parity_sign = "-" if configuration.parity == -1 else "+"

    print(f"Configuration: {configuration.label}")
    print(f"Parity: {configuration.parity_name}")
    print(f"Microstates: {microstates}")
    print("\nLS terms and fine-structure levels:")
    for term in terms:
        suffix = configuration.parity_symbol
        occurrence_text = f" x{term.occurrences}" if term.occurrences > 1 else ""
        print(
            f"  {term.label}{suffix}{occurrence_text}  "
            f"(L={term.L}, S={format_half(term.two_s)}, "
            f"states/occurrence={term.degeneracy})"
        )
        for two_j in j_values(term):
            g_value = lande_g(term, two_j)
            g_text = "n/a" if g_value is None else f"{g_value:.6g}"
            print(
                f"      J={format_half(two_j):>3}  "
                f"2J+1={two_j + 1:<2}  g_J={g_text}"
            )

    print("\nExpected J^parity level blocks:")
    for two_j, count in j_block_counts(terms).items():
        print(
            f"  J={format_half(two_j):>3}{parity_sign}  "
            f"levels={count:<3}  magnetic sublevels={count * (two_j + 1)}"
        )

    hund, hund_note = hund_ground(configuration, terms)
    if hund is not None:
        print(
            f"\nHund estimate: {hund.term.label}{configuration.parity_symbol}_"
            f"{format_half(hund.two_j)} ({hund.filling})"
        )
    elif hund_note:
        print(f"\nHund estimate: unavailable ({hund_note})")

    if stats:
        term_states = sum(term.occurrences * term.degeneracy for term in terms)
        level_states = sum(
            count * (two_j + 1) for two_j, count in j_block_counts(terms).items()
        )
        expected_microstates = math.prod(
            math.comb(shell.capacity, shell.occupancy)
            for shell in configuration.subshells
        )
        print("\nConsistency checks:")
        print(f"  determinant count:       {microstates}")
        print(f"  combinatorial count:     {expected_microstates}")
        print(f"  LS term degeneracies:    {term_states}")
        print(f"  J-level degeneracies:    {level_states}")
        print("  status:                  OK")


def print_parentage(configuration: Configuration, paths: Sequence[ParentagePath]) -> None:
    print("\nParentage (open subshells coupled left to right):")
    if not paths:
        print("  unavailable: fewer than two open subshells")
        return
    for path in paths:
        genealogy = " -> ".join(path.intermediate_terms)
        suffix = configuration.parity_symbol
        if genealogy:
            genealogy = genealogy[: -len(path.final_label)] + path.final_label + suffix
        occurrence_text = f" x{path.occurrences}" if path.occurrences > 1 else ""
        print(f"  {' x '.join(path.parents)} => {genealogy}{occurrence_text}")


def print_jj(configuration: Configuration, analysis: JjAnalysis) -> None:
    print("\nRelativistic subshell occupations and J-coupled CSFs:")
    for row in analysis.rows:
        levels = ", ".join(
            f"J={format_half(two_j)}" + (f" x{count}" if count > 1 else "")
            for two_j, count in row.levels
        )
        print(f"  {row.label}: {levels}  ({row.microstates} magnetic sublevels)")

    parity_sign = "-" if configuration.parity == -1 else "+"
    jj_summary = "  ".join(
        f"J={format_half(two_j)}{parity_sign}: {count}"
        for two_j, count in analysis.jj_census
    )
    ls_summary = "  ".join(
        f"J={format_half(two_j)}{parity_sign}: {count}"
        for two_j, count in analysis.ls_census
    )
    print(f"  Levels per J (jj): {jj_summary}")
    print(f"  Levels per J (LS): {ls_summary}")
    print("  cross-check:       OK")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Enumerate LS terms and J-parity levels as a preflight check for "
            "atomic-structure calculations. No energies are predicted."
        )
    )
    parser.add_argument(
        "configuration",
        nargs="?",
        help="configuration such as '2p2', '2p1 3s1', or shorthand 'p2'",
    )
    parser.add_argument(
        "--stats", "-s", action="store_true", help="show independent state-count checks"
    )
    parser.add_argument(
        "--transitions",
        metavar="CONFIGURATION",
        help="list direct pure-LS E1 transitions to another configuration",
    )
    parser.add_argument(
        "--json", action="store_true", help="emit machine-readable JSON"
    )
    parser.add_argument(
        "--parentage",
        action="store_true",
        help="show how terms of multiple open subshells couple",
    )
    parser.add_argument(
        "--jj",
        action="store_true",
        help="show relativistic occupations and J-coupled CSF counts",
    )
    parser.add_argument(
        "--nist-spectrum",
        metavar="SPECTRUM",
        help="include a NIST ASD levels link, for example 'Ti I'",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    config_text = arguments.configuration
    if config_text is None:
        config_text = input("Enter electron configuration (for example, 2p2): ").strip()

    try:
        configuration = parse_configuration(config_text)
        terms = terms_for_configuration(configuration)
        parentage_result = (
            parentage_paths(configuration, terms, terms_for_subshell)
            if arguments.parentage
            else ()
        )
        jj_result = jj_analysis(configuration, terms) if arguments.jj else None
        final_configuration = None
        final_terms = None
        transition_result = None
        if arguments.transitions:
            final_configuration = parse_configuration(arguments.transitions)
            final_terms = terms_for_configuration(final_configuration)
            transition_result = allowed_e1_transitions(
                configuration, terms, final_configuration, final_terms
            )
    except ConfigurationError as error:
        parser.error(str(error))

    if arguments.json:
        document = {"initial": analysis_dict(configuration, terms)}
        if arguments.parentage:
            document["initial"]["parentage"] = parentage_dict(
                parentage_result, configuration.parity
            )
        if jj_result is not None:
            document["initial"]["jj_coupling"] = jj_dict(
                jj_result, configuration.parity
            )
        if final_configuration is not None and final_terms is not None:
            document["final"] = analysis_dict(final_configuration, final_terms)
            document["e1_transitions"] = transition_dict(
                transition_result, configuration.parity, final_configuration.parity
            )
        if arguments.nist_spectrum:
            document["nist_asd_levels"] = nist_levels_url(arguments.nist_spectrum)
        print(json.dumps(document, indent=2, sort_keys=True))
        return 0

    print_analysis(configuration, terms, arguments.stats)
    if arguments.parentage:
        print_parentage(configuration, parentage_result)
    if jj_result is not None:
        print_jj(configuration, jj_result)
    if final_configuration is not None and final_terms is not None:
        print("\nTransition target:")
        print_analysis(final_configuration, final_terms, arguments.stats)
        print("\nDirect pure-LS E1 candidates:")
        if transition_result.reason:
            print(f"  none: {transition_result.reason}")
        else:
            for transition in transition_result.transitions:
                initial = (
                    f"{transition.initial_term.label}{configuration.parity_symbol}_"
                    f"{format_half(transition.initial_two_j)}"
                )
                final = (
                    f"{transition.final_term.label}{final_configuration.parity_symbol}_"
                    f"{format_half(transition.final_two_j)}"
                )
                suffix = (
                    f" x{transition.combinations} term occurrences"
                    if transition.combinations > 1
                    else ""
                )
                print(f"  {initial} -> {final}{suffix}")

    if arguments.nist_spectrum:
        print(f"\nNIST ASD levels: {nist_levels_url(arguments.nist_spectrum)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
