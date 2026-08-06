# Multiplet Generator

`multiplet_generator.py` enumerates the LS terms and fine-structure level
symmetries expected from an electron configuration. Its main use is a first
pass before GRASP, Cowan-code, or similar atomic-structure calculations:

- Which LS terms should this configuration produce?
- How many times does each term occur?
- Which J-parity blocks should be present, and how many levels belong to each?
- Do determinant, LS-term, and J-level state counts agree?
- What ground level do Hund's rules suggest for a single open subshell?

The program does not calculate term or level energies. A configuration alone
does not supply the radial integrals, configuration interaction, or
spin-orbit matrix elements needed for physical energies.

## Quick start

The script uses only the Python standard library.

```bash
python multiplet_generator.py 2p2
python multiplet_generator.py 3d3 --stats
python multiplet_generator.py "2p1 3s1" --parentage --jj
```

Example for `2p2`:

```text
Configuration: 2p2
Parity: even
Microstates: 15

LS terms and fine-structure levels:
  3P  (L=1, S=1, states/occurrence=9)
      J=  0  2J+1=1   g_J=n/a
      J=  1  2J+1=3   g_J=1.5
      J=  2  2J+1=5   g_J=1.5
  1D  (L=2, S=0, states/occurrence=5)
      J=  2  2J+1=5   g_J=1
  1S  (L=0, S=0, states/occurrence=1)
      J=  0  2J+1=1   g_J=n/a

Expected J^parity level blocks:
  J=  0+  levels=2    magnetic sublevels=2
  J=  1+  levels=1    magnetic sublevels=3
  J=  2+  levels=2    magnetic sublevels=10

Hund estimate: 3P_0 (less than half-filled)
```

`g_J=n/a` for `J=0` because the usual Landé expression has a zero
denominator; the linear magnetic moment is zero.

## Configuration syntax

Use one or more `[principal quantum number][orbital][occupancy]` tokens:

```text
2p2
2p1 3s1
3d5 4s1
```

The principal quantum number may be omitted when only angular structure
matters:

```text
p2
d3
s1p4
```

Electrons within every `nl` subshell are equivalent and antisymmetrized.
Different subshells are coupled after their individual determinant spaces are
built. Thus `2p2` and `p2` have the same LS terms, while `2p1 3p1` describes
two non-equivalent groups and has a larger term set.

Repeated tokens for the same subshell are merged before electron-hole
reduction, so `p4p1` is treated as `p5`. Occupancies above the subshell
capacity and impossible combinations such as `1p1` are rejected.

The parser accepts the spectroscopic orbital sequence
`s p d f g h i k l m n o q`; `j` is omitted by convention.

## State-count checks

Use `--stats` to show four independently assembled counts:

```bash
python multiplet_generator.py 3d3 --stats
```

The program compares:

1. The determinant weight distribution.
2. The product of binomial subshell counts.
3. The sum of `(2L+1)(2S+1)` over LS terms.
4. The sum of `2J+1` over fine-structure levels.

All four must agree. A disagreement is treated as an implementation error,
including for configurations with multiple non-equivalent subshells.

## Parentage across open subshells

Use `--parentage` when a configuration has two or more open subshells:

```bash
python multiplet_generator.py "3d2 4s1" --parentage
python multiplet_generator.py "3d3 4s1" --parentage
```

The program first derives the allowed terms of each antisymmetrized subshell,
then couples their L and S values from left to right. Each genealogy retains
the number of occurrences inherited from recurring parent terms. For example,
the two `2D` terms of `3d3` produce two occurrences of both `3D` and `1D`
when coupled to `4s1(2S)`.

This labels parent subshell terms and intermediate couplings. It does not
assign seniority or Nielson-Koster labels to recurring terms within one
equivalent-electron subshell.

## Relativistic occupations and GRASP CSFs

Use `--jj` to split every open `nl` subshell over `j=l-1/2` and `j=l+1/2`:

```bash
python multiplet_generator.py 2p2 --jj
python multiplet_generator.py 4f7 --jj
python multiplet_generator.py "2p1 3s1" --jj
```

For each allowed relativistic occupation, the program counts antisymmetrized
`M_J` determinants and obtains the number of J-coupled CSFs from the adjacent
weight difference. The summed levels-per-J census is independently compared
with the LS decomposition and must match exactly.

For `2p2`, for example, the relativistic rows are
`2p_1/2^2`, `2p_1/2^1 2p_3/2^1`, and `2p_3/2^2`. Their combined census is
two `J=0` levels, one `J=1` level, and two `J=2` levels, matching the LS side.

## Pure-LS E1 candidates

`--transitions` compares two configurations and lists direct electric-dipole
candidates:

```bash
python multiplet_generator.py 2p1 --transitions 3s1
python multiplet_generator.py "2p1 3s1" --transitions "2p2"
```

The filter requires a one-electron promotion with `Delta l = +/-1`, a parity
change, `Delta S = 0`, `Delta L = 0, +/-1` excluding `L=0 -> 0`, and
`Delta J = 0, +/-1` excluding `J=0 -> 0`.

These are pure-LS, single-configuration candidates. Configuration mixing can
give nominally forbidden lines nonzero intensity, and this program does not
calculate line strengths.

## JSON output

Use `--json` for scripts, notebooks, or comparison tooling:

```bash
python multiplet_generator.py 3d3 --json
python multiplet_generator.py 2p1 --transitions 3s1 --json
```

The JSON includes terms, occurrence counts, J levels, Landé factors, J-parity
block counts, microstate counts, and Hund guidance. `--parentage` and `--jj`
add their corresponding records to the JSON document.

## NIST ASD comparison link

If an element and ionization stage are known, the CLI can include a direct
NIST Atomic Spectra Database levels query:

```bash
python multiplet_generator.py 3d2 --nist-spectrum "Ti I"
```

This prints a query link; it does not scrape or reinterpret ASD assignments.
Automatic matching would be unreliable for mixed configurations and levels
whose names represent only the leading eigenvector component.

NIST references used for the conventions in this program:

- [Atomic states, shells, configurations, and parity](https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy-10)
- [Allowed terms and recurring-term labels](https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy-1)
- [LS-coupling transition selection rules](https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy)
- [NIST ASD energy-level search help](https://physics.nist.gov/PhysRefData/ASD/Html/levelshelp.html)

## Relationship to GRASP output

The default `J^parity` summary gives the physical level count implied by the
nonrelativistic configuration. The `--jj` rows give the CSF count for each
specific relativistic occupation and reproduce the same total J census when
all occupations belonging to the nonrelativistic configuration are included.

Several distinctions matter:

- A per-occupation GRASP CSF count is a basis dimension. The summed census
  over every relativistic split of one nonrelativistic configuration also
  equals the number of levels supplied by the LS decomposition.
- GRASP uses relativistic subshells and jj-coupled CSFs, while this program
  classifies the nonrelativistic configuration in pure LS coupling.
- Configuration interaction mixes levels with the same J and parity.
- Strong relativistic mixing can make an LS label only an approximate name.
- Repeated terms such as the two `2D` terms of `d3` are counted separately.
  `--parentage` carries that count through multi-subshell coupling, while
  seniority labels within `d3` remain outside the current scope.

The counts are therefore a symmetry and bookkeeping check. Energy ordering,
mixing coefficients, and transition rates belong to the atomic-structure
calculation.

## Testing

Run the standard-library test suite with:

```bash
python -m unittest -v
```

The suite pins known `p2`, `d2`, and `d3` term sets; explicit-subshell and
parser behavior; Hund ground levels; E1 filters; JSON output; parentage; a
brute-force differential check of the dynamic j-subshell counter; and
exhaustive LS, jj, state-count, and electron-hole checks for every occupancy
from `s` through `f`.

## Removed energy options

Earlier development versions accepted `--so`, `--zeta`, `--rs_scale`, `--a`,
and `--b`. Those options used a generic `S(S+1)`/`L(L+1)` expression that can
even violate Hund ordering, so they have been removed rather than deprecated.
Another experimental rewrite added limited Racah-energy tables and
octahedral-field reductions. Those remain in Git history but were not folded
into the preflight CLI because their physical scope differs from a free-ion
symmetry census.
Use GRASP, Cowan code, Quanty, or another atomic-structure program when actual
energies or mixing are required.
