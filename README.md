# Multiplet Generator

A Python-based tool to generate LS multiplets from an electron configuration using the LS coupling (Russell–Saunders coupling) scheme. This program calculates the full microstate distribution for a given configuration, performs an iterative "peeling" algorithm to extract the LS multiplets, and optionally provides a consistency check against the total number of unperturbed microstates.

## Features

- LS Coupling Calculation:
  - Computes possible LS multiplets (term symbols such as 3P, 1D, 1S, etc.) from the specified electron configuration.

- Equivalent & Non-Equivalent Electrons:
  - Treats tokens without a specified principal quantum number as equivalent electrons (with antisymmetrization) and tokens with a number (e.g., 1s) as non-equivalent.

- Overall Parity:
  - Computes the overall parity of the configuration using the rule (−1)∑(𝑙×count)

- Consistency Check (Statistics Mode):
  - When the --stats or -s flag is provided, the tool displays the total number of microstates (i.e. unperturbed states) and the sum of the degeneracies of the LS multiplets, so you can verify that the LS decomposition accounts for all microstates.

## Usage

Run the script from the command line with the desired configuration string:

```bash
python multiplet_generator.py s1p4
```

For additional statistics and consistency checks, add the --stats flag:

```bash
python multiplet_generator.py s1s1 --stats
```

Display help message:

```bash
python multiplet_generator.py --help
```

## Configuration Format

The configuration string is composed of tokens in the form:

```
[principal quantum number][orbital letter][electron count]
```

- Without Principal Quantum Number:
  - Tokens like p2 denote 2 equivalent electrons in a p-shell. Equivalent electrons are antisymmetrized.

- With Principal Quantum Number:
  - Tokens like 1s1 denote an electron in a unique shell (non-equivalent). Multiple non-equivalent tokens (e.g., s1s1) are not merged.

- Combined Example:
  - "1sp2" indicates one electron in the 1s orbital and 2 electrons in a p orbital.

## How It Works

1. Parsing & Tokenization:
The code parses the input string into tokens. Tokens with an explicit principal quantum number are treated as non-equivalent, while tokens without one are merged.

2. Microstate Distribution:
For each shell, a microstate distribution is computed:

- Equivalent electrons: Use Slater determinants to ensure antisymmetry.
- Non-equivalent electrons: Use independent single-electron state counting.

3. Distribution Combination:
The individual distributions from each shell are combined via convolution to obtain the overall microstate distribution.

4. Multiplet Peeling:
A peeling algorithm iteratively subtracts multiplet blocks (each representing a full LS multiplet with degeneracy (2L+1)(2S+1)) from the overall distribution. This produces the LS multiplets along with the number of times each multiplet is subtracted.

5. Parity Calculation:
The overall configuration parity is computed using (−1) (sum of count * l)
  and a superscript ° is appended if the parity is odd.

6. Output & Consistency Check:

The tool prints:

- Overall parity.
- A list of LS multiplets with their term symbols, counts (from peeling), and expected degeneracies.
- (Optional) A consistency check showing the total number of microstates versus the sum of LS multiplet degeneracies.

## Requirements
Python 3 (no additional libraries required)

## Future Enhancements

- Adding a verbose/debug mode to output intermediate microstate distributions.
- Graphical output of the microstate distribution (e.g., a histogram of ML vs MS)
- Extending the algorithm to support alternative coupling schemes (such as j–j coupling).



