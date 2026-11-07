# Reproducibility — The L’vov–Kaplansky Conjecture for Standard Polynomials

**Y. T. Miah · 7 November 2026**

The article is built from `main.tex`, `odd_proof.tex`, `even_proof.tex`, `verification.tex` and `references.bib`. The accompanying programs check specified matrix-unit identities, small exceptional orders, and exact integer witnesses for admissible cores. The general proofs and external generic-nullity dependencies are given in the article.

## Requirements and odd-standard checks

Use Python 3.10 or later with the packages in `requirements.txt` (`sympy==1.14.0` and NumPy). From the manuscript directory run:

```bash
python -m pip install -r requirements.txt
python code/verify_standard.py --output standard_rerun.json
python code/verify_standard.py --extended --output standard_extended_rerun.json
```

`code/verify_standard.py` calculates signed alternators by a full-matrix subset recurrence and compares selected results with the cyclic-homogeneous recurrence in `code/graded.py` and Euler-trail enumeration. The tests cover local paths and cycles, periodic parameters, exceptional finite orders and sample Jordan constructions. The examples range over matrix orders 3–8, but direct full-alternator evaluation is checked only in orders 3–5 (26 cases); the 30 cases in orders 6–8 retain `direct_full_evaluation: false` and record construction parameters rather than the full input matrices. The default run checks 371 local path instances, 740 cycle instances and 18 periodic pairs; `--extended` checks 636, 3,190 and 25, respectively.

Compare the **parsed JSON values** of `standard_rerun.json` and `standard_extended_rerun.json` with `certificates/standard_certificates.json` and `certificates/standard_certificates_extended.json`. **Do not use `python -O`**: assertions are part of the tests.

## Uniform even-degree seed

The program `code/verify_minimal_even.py` checks the signed matrix-unit formulas and computes the full differential in selected finite orders:

```bash
python code/verify_minimal_even.py --max-n 8 --rank-max-n 6 --output minimal_even_rerun.json
```

The symbolic proof for all orders is in the article; the program checks the finite sizes selected by its command-line arguments.

## Seven bordered-core certificates

The certificate files are:

```text
certificates/border_r4_N5.json
certificates/border_r4_N6.json
certificates/border_r4_N7.json
certificates/border_r5_N6.json
certificates/border_r5_N7.json
certificates/border_r5_N8.json
certificates/border_r5_N9.json
```

Each record contains an explicit integer tuple and column, an endpoint determinant and a selected full-size normal-map minor. The normal map has `N*N` columns; a nonzero `N*N` minor proves injectivity. The saved-input checker computes the endpoint determinant exactly and independently reconstructs all normal-map columns. It computes the selected full-size minor **modulo 101 and 103**, comparing both residues with those of the recorded integer minor; a nonzero residue certifies that this minor is nonzero in characteristic zero. It does not independently recompute the entire integer minor. The separate constructor computes the integer determinant for a regenerated case. The uniform seed at core size `N=r`, the seven finite certificates and the additive-core theorem establish all positive-size cases in degrees eight and ten; the zero cases follow from the Amitsur–Levitzki identity.

To verify the saved inputs, run the following for **each** of the seven filenames above, changing both the input filename and output name:

```bash
python code/verify_saved_border.py \
  --input certificates/border_r4_N5.json \
  --output checked_border_r4_N5.json
```

`code/verify_saved_border.py` reconstructs all normal-map columns by a first-letter recurrence and checks the endpoint, minor residues and bordered rank-one square-zero value. The constructor `code/verify_border_even.py` uses a different subset recurrence and can regenerate a saved case from its deterministic search settings, for example:

```bash
python code/verify_border_even.py --case 4 5 --output regenerated_r4_N5.json
```

Constructing new witnesses can take considerably longer than verifying the seven saved inputs. The two programs use different recurrences, allowing the saved inputs to be checked separately from witness construction. The constructor accepts only the seven documented `(r,N)` cases.

## Building the article

From this directory, with a LaTeX installation containing `amsart`, `amsmath`, `amssymb`, `mathtools`, `microtype`, `hyperref`, `enumitem`, `booktabs` and `lmodern`, run:

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex8 main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The root `main.tex` includes the two proof files and the verification appendix. The resulting `main.pdf` is the manuscript PDF.
