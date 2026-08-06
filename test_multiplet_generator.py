"""Regression tests for term enumeration and the command-line interface.

The cases pin known LS terms, state-count invariants, and parser boundaries.
"""

import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

import multiplet_generator as mg


SCRIPT = Path(__file__).with_name("multiplet_generator.py")


def term_counts(configuration_text):
    configuration = mg.parse_configuration(configuration_text)
    return {
        term.label: term.occurrences
        for term in mg.terms_for_configuration(configuration)
    }


class ConfigurationTests(unittest.TestCase):
    def test_principal_quantum_number_does_not_change_equivalence(self):
        self.assertEqual(term_counts("p2"), term_counts("2p2"))
        self.assertEqual(term_counts("2p2"), {"3P": 1, "1D": 1, "1S": 1})

    def test_repeated_subshells_merge_before_hole_reduction(self):
        configuration = mg.parse_configuration("p4p1")
        self.assertEqual(configuration.subshells[0].occupancy, 5)
        self.assertEqual(term_counts("p4p1"), {"2P": 1})

    def test_multiple_principal_subshells_remain_distinct(self):
        self.assertEqual(
            term_counts("2p1 3p1"),
            {"3D": 1, "3P": 1, "3S": 1, "1D": 1, "1P": 1, "1S": 1},
        )

    def test_invalid_configurations_are_rejected(self):
        for text in ("", "nonsense", "1sp2", "1p1", "2p7", "j1", "p0"):
            with self.subTest(text=text):
                with self.assertRaises(mg.ConfigurationError):
                    mg.parse_configuration(text)

    def test_extended_orbital_letters_are_consistent(self):
        configuration = mg.parse_configuration("l1")
        self.assertEqual(configuration.subshells[0].l, 8)
        self.assertEqual(term_counts("l1"), {"2L": 1})


class TermEnumerationTests(unittest.TestCase):
    def test_known_equivalent_electron_terms(self):
        self.assertEqual(term_counts("p2"), {"3P": 1, "1D": 1, "1S": 1})
        self.assertEqual(
            term_counts("d2"),
            {"3F": 1, "3P": 1, "1G": 1, "1D": 1, "1S": 1},
        )
        self.assertEqual(
            term_counts("d3"),
            {"4F": 1, "4P": 1, "2H": 1, "2G": 1, "2F": 1, "2D": 2, "2P": 1},
        )

    def test_all_s_through_f_occupancies_close_exactly(self):
        for l_value in range(4):
            capacity = 2 * (2 * l_value + 1)
            for occupancy in range(capacity + 1):
                with self.subTest(l=l_value, occupancy=occupancy):
                    distribution = mg.shell_microstate_distribution(l_value, occupancy)
                    terms = mg.decompose_ls(distribution)
                    states = sum(
                        term.occurrences * term.degeneracy for term in terms
                    )
                    self.assertEqual(states, math.comb(capacity, occupancy))

    def test_electron_hole_term_equivalence(self):
        for l_value in range(4):
            capacity = 2 * (2 * l_value + 1)
            for occupancy in range(capacity // 2 + 1):
                with self.subTest(l=l_value, occupancy=occupancy):
                    electron_terms = mg.decompose_ls(
                        mg.shell_microstate_distribution(l_value, occupancy)
                    )
                    hole_terms = mg.decompose_ls(
                        mg.shell_microstate_distribution(l_value, capacity - occupancy)
                    )
                    self.assertEqual(electron_terms, hole_terms)

    def test_non_equivalent_p_s_configuration(self):
        configuration = mg.parse_configuration("2p1 3s1")
        self.assertEqual(configuration.parity_name, "odd")
        self.assertEqual(term_counts("2p1 3s1"), {"3P": 1, "1P": 1})


class LevelTests(unittest.TestCase):
    def test_triplet_p_levels_and_lande_factors(self):
        triplet_p = mg.Term(L=1, two_s=2, occurrences=1)
        self.assertEqual(mg.j_values(triplet_p), (0, 2, 4))
        self.assertIsNone(mg.lande_g(triplet_p, 0))
        self.assertAlmostEqual(mg.lande_g(triplet_p, 2), 1.5)
        self.assertAlmostEqual(mg.lande_g(triplet_p, 4), 1.5)

    def test_p2_j_block_counts(self):
        configuration = mg.parse_configuration("p2")
        blocks = mg.j_block_counts(mg.terms_for_configuration(configuration))
        self.assertEqual(blocks, {0: 2, 2: 1, 4: 2})

    def test_hund_ground_levels(self):
        expected = {
            "p2": ("3P", 0),
            "p4": ("3P", 4),
            "d2": ("3F", 4),
            "d5": ("6S", 5),
            "d8": ("3F", 8),
        }
        for text, (label, two_j) in expected.items():
            with self.subTest(text=text):
                configuration = mg.parse_configuration(text)
                terms = mg.terms_for_configuration(configuration)
                hund, note = mg.hund_ground(configuration, terms)
                self.assertIsNone(note)
                self.assertEqual((hund.term.label, hund.two_j), (label, two_j))

    def test_hund_guidance_stops_at_multiple_open_subshells(self):
        configuration = mg.parse_configuration("2p1 3s1")
        terms = mg.terms_for_configuration(configuration)
        hund, note = mg.hund_ground(configuration, terms)
        self.assertIsNone(hund)
        self.assertIn("one open subshell", note)


class TransitionTests(unittest.TestCase):
    def test_direct_p_to_s_e1_candidates(self):
        initial = mg.parse_configuration("2p1")
        final = mg.parse_configuration("3s1")
        result = mg.allowed_e1_transitions(
            initial,
            mg.terms_for_configuration(initial),
            final,
            mg.terms_for_configuration(final),
        )
        pairs = {
            (transition.initial_two_j, transition.final_two_j)
            for transition in result.transitions
        }
        self.assertIsNone(result.reason)
        self.assertEqual(pairs, {(1, 1), (3, 1)})

    def test_e1_rejects_same_parity(self):
        initial = mg.parse_configuration("2p1")
        final = mg.parse_configuration("3p1")
        result = mg.allowed_e1_transitions(
            initial,
            mg.terms_for_configuration(initial),
            final,
            mg.terms_for_configuration(final),
        )
        self.assertEqual(result.transitions, ())
        self.assertIn("parity", result.reason)


class CommandLineTests(unittest.TestCase):
    def run_cli(self, *arguments):
        return subprocess.run(
            [sys.executable, str(SCRIPT), *arguments],
            capture_output=True,
            check=False,
            text=True,
        )

    def test_json_output(self):
        completed = self.run_cli("2p2", "--json", "--nist-spectrum", "Ti I")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        document = json.loads(completed.stdout)
        self.assertEqual(document["initial"]["microstates"], 15)
        self.assertEqual(document["initial"]["hund_ground"]["term"], "3P")
        self.assertIn("spectrum=Ti+I", document["nist_asd_levels"])

    def test_invalid_input_has_nonzero_exit(self):
        completed = self.run_cli("2p7")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("exceeds the subshell capacity", completed.stderr)

    def test_removed_energy_option_is_rejected(self):
        completed = self.run_cli("p2", "--rs_scale=1000")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("unrecognized arguments", completed.stderr)


if __name__ == "__main__":
    unittest.main()
