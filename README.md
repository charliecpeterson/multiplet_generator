# Multiplet Generator

LS-coupling term lookup for electron configurations, so I don't have to keep
finding Table 13-2 in Slater's *Quantum Theory of Atomic Structure*. Also a
sanity checker for GRASP / DIRAC runs: the `--jj` output reproduces the
relativistic-subshell CSF bookkeeping (levels per J and parity) that
`rcsfgenerate`/`rlevels` report, cross-checked against the LS decomposition.
Single file, stdlib only.

For a configuration it prints the allowed LS terms with J levels, term
parentage for multi-shell configurations, the Hund's-rules ground state with
Landé g and effective magnetic moments, and on request: the M_L/M_S
microstate table, jj-coupling splits, interval-rule fine structure,
octahedral-field term splitting, and the classic diagonal-sum term energies
for p²/p³/d²/d³ (and hole partners).

## Usage

```bash
python multiplet_generator.py d3
python multiplet_generator.py 2p1 3p1            # two nonequivalent p electrons
python multiplet_generator.py 1s2.2s2.2p3       # nitrogen; . , / or space separate tokens
python multiplet_generator.py 3d2 4s1           # multi-shell: terms come with parentage
python multiplet_generator.py p2 --grid         # M_L/M_S table (Slater Fig. 13-3 style)
python multiplet_generator.py f7 --jj           # relativistic CSF splits, levels per J
python multiplet_generator.py d3 --zeta 273     # fine structure via interval rule
python multiplet_generator.py d2 --oh           # weak-field splitting in Oh
python multiplet_generator.py d2 --energies     # Racah expressions, symbolic
python multiplet_generator.py d3 --B 918 --C 3850   # evaluated, in the unit of B and C
python multiplet_generator.py p2 --F2 1500      # p^n uses Slater-Condon F2 instead
```

## Configuration format

Tokens are `[n]<letter><count>`: `p2`, `3d4`, `2p1`. Rules:

- Electrons within a token share a subshell and are equivalent
  (antisymmetrized, Pauli applies).
- Different tokens are independent subshells, coupled by convolution of
  their microstate distributions. So `2p1 3p1` is the nonequivalent-pair
  case, and `p1 p1` (no n given) also means two distinct p subshells.
- Tokens with the same explicit n and letter merge: `2p1 2p1` ≡ `2p2`.
- Concatenated multi-token strings like `2p13p1` are ambiguous and
  rejected; separate tokens with space, `.`, `,`, or `/`.

## What it computes

1. **Terms.** Slater determinants per subshell tallied by (M_L, M_S),
   convolved across subshells; the number of (L, S) terms is
   `n(L,S) − n(L+1,S) − n(L,S+1) + n(L+1,S+1)` — Slater's Sec. 13-2
   bookkeeping in closed form.
2. **Parentage** (automatic for ≥2 open subshells). Each subshell's own
   terms are coupled left to right with the triangle rule, giving the
   `3d²(³F)4s → ⁴F, ²F` genealogy used by NIST and GRASP labels. The
   coupled totals are asserted against the microstate count.
3. **Hund's rules** (single open subshell): ground term, ground J,
   normal/inverted multiplet, Landé g per level, and μ_eff both spin-only
   (typical quenched 3d complex) and g_J√(J(J+1)) (free ion, 4f/5f).
4. **jj coupling** (`--jj`). Splits each open nl shell over j = l±1/2,
   enumerates the relativistic configurations, and peels J values from the
   M_J distributions. The levels-per-J census must match the LS side, and
   both are printed — this is the GRASP sanity check.
5. **Fine structure** (`--zeta Z`). Landé interval rule on the ground term,
   λ = ±ζ/2S, sign from the filling. First-order only: real f-ion levels
   deviate (Pr³⁺ ³H₅ comes out ~1850 vs. ~2150 cm⁻¹ observed).
6. **Oh splitting** (`--oh`). Reduction of χ_L onto the irreps of O with
   g/u from the configuration parity — the weak-field limit of
   Tanabe–Sugano (³F → ³A₂g + ³T₂g + ³T₁g and so on).
7. **Hole equivalence.** l^n and l^(4l+2−n) have identical terms and
   electrostatic splittings; noted, and used for energy lookups (d⁸ → d²).

## Term energies

`--energies` prints the diagonal-sum-rule results (Slater Secs. 13-1–13-3):
p² and p³ in Condon–Shortley F2 (= F²/25), d² and d³ in Racah A, B, C.
The contribution common to all terms of the configuration (F0 terms,
multiples of A) is omitted since it cancels in splittings. The d³ ²D pair
comes from the 2×2 secular problem, hence the square root; ²P and ²H of d³
are degenerate at this level of treatment.

Coefficients are the standard tabulated ones (Condon & Shortley; Griffith,
*The Theory of Transition-Metal Ions*, Table 4.6). Cross-checks: term lists
for p^n, d^n, f² reproduce Slater's Table 13-2 including repeated-term
counts (two ²D in d³, three ²D in d⁵); f⁷ gives the known 119 terms and
327 levels, with the jj and LS per-J censuses agreeing; the f² ³H₄ Landé g
is 0.800 and μ_eff 3.58 μ_B (Pr³⁺ textbook values).

No B/C values are built in — free-ion Racah parameters vary by source, so
pass your own (any unit; output splittings come out in the same unit).

## Possible extensions

- Selection-rule checker between two configurations (Laporte, ΔS, ΔL, ΔJ).
- Generic diagonal-sum energies from c^k coefficients for terms that appear
  once, instead of hardcoded tables.
- Restricted jj occupations (active-space style) to mirror a specific
  rcsfgenerate input rather than the full expansion.
