"""GRASP-facing parentage, jj-coupling, and transition analyses.

The routines operate on the small configuration and term objects supplied by
the CLI module while keeping relativistic bookkeeping in one place.
"""

import itertools
import math
from dataclasses import dataclass
from typing import Callable, Dict, Mapping, Optional, Protocol, Sequence, Tuple


class ShellLike(Protocol):
    n: Optional[int]
    l: int
    occupancy: int
    label: str
    name: str
    is_open: bool


class ConfigurationLike(Protocol):
    subshells: Tuple[ShellLike, ...]
    electron_count: int
    parity: int


class TermLike(Protocol):
    L: int
    two_s: int
    occurrences: int
    label: str


TermProvider = Callable[[ShellLike], Sequence[TermLike]]


@dataclass(frozen=True)
class ParentagePath:
    parents: Tuple[str, ...]
    intermediate_terms: Tuple[str, ...]
    final_L: int
    final_two_s: int
    final_label: str
    occurrences: int


@dataclass(frozen=True)
class RelativisticOccupation:
    label: str
    two_j: int
    occupancy: int


@dataclass(frozen=True)
class JjRow:
    label: str
    occupations: Tuple[RelativisticOccupation, ...]
    levels: Tuple[Tuple[int, int], ...]
    microstates: int


@dataclass(frozen=True)
class JjAnalysis:
    rows: Tuple[JjRow, ...]
    jj_census: Tuple[Tuple[int, int], ...]
    ls_census: Tuple[Tuple[int, int], ...]

    @property
    def consistent(self) -> bool:
        return self.jj_census == self.ls_census


@dataclass(frozen=True)
class Transition:
    initial_term: TermLike
    initial_two_j: int
    final_term: TermLike
    final_two_j: int

    @property
    def combinations(self) -> int:
        return self.initial_term.occurrences * self.final_term.occurrences


@dataclass(frozen=True)
class TransitionResult:
    transitions: Tuple[Transition, ...]
    reason: Optional[str]


@dataclass(frozen=True)
class _ParentageState:
    L: int
    two_s: int
    parents: Tuple[str, ...]
    intermediate_terms: Tuple[str, ...]
    occurrences: int


def _format_half(two_value: int) -> str:
    if two_value % 2 == 0:
        return str(two_value // 2)
    return f"{two_value}/2"


def _coupled_terms(L1: int, two_s1: int, L2: int, two_s2: int):
    for L in range(abs(L1 - L2), L1 + L2 + 1):
        for two_s in range(abs(two_s1 - two_s2), two_s1 + two_s2 + 1, 2):
            yield L, two_s


def parentage_paths(
    configuration: ConfigurationLike,
    terms: Sequence[TermLike],
    term_provider: TermProvider,
) -> Tuple[ParentagePath, ...]:
    """Couple open-subshell terms left to right and retain each genealogy."""
    open_shells = [shell for shell in configuration.subshells if shell.is_open]
    if len(open_shells) < 2:
        return ()

    first_shell = open_shells[0]
    states = [
        _ParentageState(
            L=term.L,
            two_s=term.two_s,
            parents=(f"{first_shell.label}({term.label})",),
            intermediate_terms=(),
            occurrences=term.occurrences,
        )
        for term in term_provider(first_shell)
    ]

    for shell in open_shells[1:]:
        updated = []
        for state in states:
            for shell_term in term_provider(shell):
                parent = f"{shell.label}({shell_term.label})"
                for L, two_s in _coupled_terms(
                    state.L, state.two_s, shell_term.L, shell_term.two_s
                ):
                    final_label = _term_label(L, two_s)
                    updated.append(
                        _ParentageState(
                            L=L,
                            two_s=two_s,
                            parents=state.parents + (parent,),
                            intermediate_terms=state.intermediate_terms
                            + (final_label,),
                            occurrences=state.occurrences
                            * shell_term.occurrences,
                        )
                    )
        states = updated

    paths = tuple(
        ParentagePath(
            parents=state.parents,
            intermediate_terms=state.intermediate_terms,
            final_L=state.L,
            final_two_s=state.two_s,
            final_label=_term_label(state.L, state.two_s),
            occurrences=state.occurrences,
        )
        for state in states
    )

    tally: Dict[Tuple[int, int], int] = {}
    for path in paths:
        key = (path.final_L, path.final_two_s)
        tally[key] = tally.get(key, 0) + path.occurrences
    expected = {(term.L, term.two_s): term.occurrences for term in terms}
    if tally != expected:
        raise RuntimeError("parentage coupling disagrees with LS decomposition")
    return paths


def _term_label(L: int, two_s: int) -> str:
    letters = "SPDFGHIKLMNOQRTUVWXYZ"
    orbital = letters[L] if L < len(letters) else f"L({L})"
    return f"{two_s + 1}{orbital}"


def j_shell_m_distribution(two_j: int, occupancy: int) -> Dict[int, int]:
    """Return the 2M_J weights for equivalent electrons in a j subshell."""
    capacity = two_j + 1
    if not 0 <= occupancy <= capacity:
        raise ValueError(f"occupancy must be between 0 and {capacity}")

    effective_occupancy = min(occupancy, capacity - occupancy)
    states: Dict[Tuple[int, int], int] = {(0, 0): 1}
    for two_m_j in range(-two_j, two_j + 1, 2):
        updated = dict(states)
        for (chosen, total_two_m_j), count in states.items():
            if chosen == effective_occupancy:
                continue
            key = (chosen + 1, total_two_m_j + two_m_j)
            updated[key] = updated.get(key, 0) + count
        states = updated

    return {
        total_two_m_j: count
        for (chosen, total_two_m_j), count in states.items()
        if chosen == effective_occupancy
    }


def _combine_m_distributions(
    distributions: Sequence[Mapping[int, int]],
) -> Dict[int, int]:
    combined = {0: 1}
    for distribution in distributions:
        updated: Dict[int, int] = {}
        for left_m, left_count in combined.items():
            for right_m, right_count in distribution.items():
                total_m = left_m + right_m
                updated[total_m] = updated.get(total_m, 0) + left_count * right_count
        combined = updated
    return combined


def extract_j_levels(distribution: Mapping[int, int]) -> Dict[int, int]:
    """Obtain J multiplicities from a symmetric 2M_J weight distribution."""
    if not distribution:
        raise ValueError("M_J distribution is empty")
    maximum_two_j = max(abs(two_m_j) for two_m_j in distribution)
    parity = next(iter(distribution)) % 2
    levels = {}
    for two_j in range(parity, maximum_two_j + 1, 2):
        count = distribution.get(two_j, 0) - distribution.get(two_j + 2, 0)
        if count < 0:
            raise RuntimeError(f"invalid jj decomposition at 2J={two_j}: {count}")
        if count:
            levels[two_j] = count
    return levels


def _jj_shell_options(shell: ShellLike):
    lower_two_j = 2 * shell.l - 1
    upper_two_j = 2 * shell.l + 1
    lower_capacity = max(0, 2 * shell.l)
    upper_capacity = 2 * shell.l + 2
    minimum_lower = max(0, shell.occupancy - upper_capacity)
    maximum_lower = min(shell.occupancy, lower_capacity)

    options = []
    for lower_occupancy in range(maximum_lower, minimum_lower - 1, -1):
        upper_occupancy = shell.occupancy - lower_occupancy
        occupations = []
        if lower_occupancy:
            occupations.append(
                RelativisticOccupation(
                    label=(
                        f"{shell.name}_{_format_half(lower_two_j)}^"
                        f"{lower_occupancy}"
                    ),
                    two_j=lower_two_j,
                    occupancy=lower_occupancy,
                )
            )
        if upper_occupancy:
            occupations.append(
                RelativisticOccupation(
                    label=(
                        f"{shell.name}_{_format_half(upper_two_j)}^"
                        f"{upper_occupancy}"
                    ),
                    two_j=upper_two_j,
                    occupancy=upper_occupancy,
                )
            )
        options.append(tuple(occupations))
    return tuple(options)


def jj_analysis(
    configuration: ConfigurationLike, terms: Sequence[TermLike]
) -> JjAnalysis:
    """Enumerate relativistic occupations and their J-coupled CSF counts."""
    open_shells = [shell for shell in configuration.subshells if shell.is_open]
    if not open_shells:
        census = ((0, 1),)
        return JjAnalysis(
            rows=(JjRow("closed shells", (), census, 1),),
            jj_census=census,
            ls_census=census,
        )

    per_shell_options = [_jj_shell_options(shell) for shell in open_shells]
    rows = []
    jj_census: Dict[int, int] = {}
    for option_set in itertools.product(*per_shell_options):
        occupations = tuple(
            occupation for shell_option in option_set for occupation in shell_option
        )
        distributions = [
            j_shell_m_distribution(occupation.two_j, occupation.occupancy)
            for occupation in occupations
        ]
        levels = extract_j_levels(_combine_m_distributions(distributions))
        level_states = sum(count * (two_j + 1) for two_j, count in levels.items())
        expected_states = math.prod(
            math.comb(occupation.two_j + 1, occupation.occupancy)
            for occupation in occupations
        )
        if level_states != expected_states:
            raise RuntimeError(
                f"jj row accounts for {level_states} of {expected_states} states"
            )
        for two_j, count in levels.items():
            jj_census[two_j] = jj_census.get(two_j, 0) + count
        rows.append(
            JjRow(
                label=" ".join(occupation.label for occupation in occupations),
                occupations=occupations,
                levels=tuple(sorted(levels.items())),
                microstates=expected_states,
            )
        )

    ls_census: Dict[int, int] = {}
    for term in terms:
        lower = abs(2 * term.L - term.two_s)
        upper = 2 * term.L + term.two_s
        for two_j in range(lower, upper + 1, 2):
            ls_census[two_j] = ls_census.get(two_j, 0) + term.occurrences

    analysis = JjAnalysis(
        rows=tuple(rows),
        jj_census=tuple(sorted(jj_census.items())),
        ls_census=tuple(sorted(ls_census.items())),
    )
    if not analysis.consistent:
        raise RuntimeError("jj and LS level censuses disagree")
    return analysis


def _direct_e1_jump(
    initial: ConfigurationLike, final: ConfigurationLike
) -> Tuple[bool, str]:
    if initial.electron_count != final.electron_count:
        return False, "the configurations have different electron counts"

    initial_counts = {(shell.n, shell.l): shell.occupancy for shell in initial.subshells}
    final_counts = {(shell.n, shell.l): shell.occupancy for shell in final.subshells}
    keys = set(initial_counts) | set(final_counts)
    changes = {
        key: final_counts.get(key, 0) - initial_counts.get(key, 0) for key in keys
    }
    removed = [key for key, change in changes.items() if change == -1]
    added = [key for key, change in changes.items() if change == 1]
    other = [change for change in changes.values() if change not in (-1, 0, 1)]
    if len(removed) != 1 or len(added) != 1 or other:
        return False, "the configurations do not differ by one electron promotion"
    if abs(removed[0][1] - added[0][1]) != 1:
        return False, "the promoted electron does not satisfy Δl = ±1"
    return True, ""


def allowed_e1_transitions(
    initial_configuration: ConfigurationLike,
    initial_terms: Sequence[TermLike],
    final_configuration: ConfigurationLike,
    final_terms: Sequence[TermLike],
) -> TransitionResult:
    if initial_configuration.parity == final_configuration.parity:
        return TransitionResult((), "E1 transitions require a parity change")

    direct_jump, reason = _direct_e1_jump(initial_configuration, final_configuration)
    if not direct_jump:
        return TransitionResult((), reason)

    transitions = []
    for initial_term in initial_terms:
        for final_term in final_terms:
            if initial_term.two_s != final_term.two_s:
                continue
            if abs(initial_term.L - final_term.L) > 1:
                continue
            if initial_term.L == 0 and final_term.L == 0:
                continue
            initial_js = range(
                abs(2 * initial_term.L - initial_term.two_s),
                2 * initial_term.L + initial_term.two_s + 1,
                2,
            )
            final_js = range(
                abs(2 * final_term.L - final_term.two_s),
                2 * final_term.L + final_term.two_s + 1,
                2,
            )
            for initial_two_j in initial_js:
                for final_two_j in final_js:
                    if abs(initial_two_j - final_two_j) > 2:
                        continue
                    if initial_two_j == 0 and final_two_j == 0:
                        continue
                    transitions.append(
                        Transition(
                            initial_term=initial_term,
                            initial_two_j=initial_two_j,
                            final_term=final_term,
                            final_two_j=final_two_j,
                        )
                    )
    return TransitionResult(tuple(transitions), None)
