# The L’vov–Kaplansky Conjecture for Standard Polynomials

**Y. T. Miah** · 7 November 2026

[Read the preprint (PDF)](main.pdf) · [LaTeX source](main.tex) · [Reproducibility](REPRODUCIBILITY.md)

**Zenodo DOI:** `10.5281/zenodo.24031055` (reserved for a draft deposit; not yet published).

## About

This preprint determines the sets of *actual values* of the unnormalized standard polynomials on full matrix algebras over algebraically closed fields of characteristic zero, in every degree and matrix order. For odd standard polynomials the image is either zero or the full matrix algebra; for even standard polynomials it is either zero or the space of traceless matrices, with the precise thresholds given in the main theorem.

The odd-degree argument uses explicit linear slices and Jordan forms. The even-degree argument develops nilpotent constructions and additive transfer, building on Qiyue Tang's work and using generic-nullity results of Brassil–Reichstein and Lin–Tang–Zhao. The manuscript gives the detailed proofs and historical attribution.

## Contents

- [`main.pdf`](main.pdf) — the 30-page manuscript.
- [`main.tex`](main.tex), [`odd_proof.tex`](odd_proof.tex), [`even_proof.tex`](even_proof.tex), [`verification.tex`](verification.tex), and [`references.bib`](references.bib) — the complete manuscript sources and bibliography.
- [`code/`](code/) — the five supplied Python verification modules and programs.
- [`certificates/`](certificates/) — nine saved JSON certificates: two standard-polynomial check records and seven even-degree bordered-core records.
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) and [`requirements.txt`](requirements.txt) — precise finite-check procedures, dependencies and manuscript build instructions.

## Building the manuscript

With a TeX installation providing the packages named in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md), run from the repository root:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex8 main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Finite verification

The provided programs check specified matrix-unit identities, selected small orders and explicit integer witnesses. The two standard-polynomial JSON records can be reproduced by the default and extended runs of `code/verify_standard.py`; the seven even-core records can be checked with `code/verify_saved_border.py`. See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for commands, coverage and the distinction between saved-input checking and witness regeneration.

## Citation

Y. T. Miah, *The L’vov–Kaplansky Conjecture for Standard Polynomials*, manuscript (2026). Reserved Zenodo DOI: `10.5281/zenodo.24031055`. Machine-readable citation data appear in [`CITATION.cff`](CITATION.cff).

## Licensing

The author's original scholarly materials and certificate data are licensed under **CC BY 4.0**; the author's original Python code is licensed under **BSD 3-Clause**. The exact file scope and full licence texts are in [`LICENSE.md`](LICENSE.md) and [`LICENSES/`](LICENSES/).
