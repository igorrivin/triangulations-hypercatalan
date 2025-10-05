# Hypercatalan Dissections of Midpoint Polygons

This repo contains:

- A **detailed LaTeX writeup** proving that dissections of the midpoint polygon
  under the constraint that every drawn segment is of the form \(M_{i,j}M_{i,k}\)
  are counted by the **order‑2 hypercatalan number**
  \[
    H(n) \;=\; \sum_{T} \;\prod_{i=1}^{n} C_{\deg_T(i)},
  \]
  where the sum is over triangulations \(T\) of a convex \(n\)-gon, \(\deg_T(i)\)
  is the number of diagonals of \(T\) incident to vertex \(i\), and \(C_k\) is the \(k\)-th Catalan number.

- A **Python program** implementing a fast dynamic program that computes \(H(n)\) exactly.

## Quick start

### 1) Compute values

```bash
python3 src/hypercatalan.py --n 10
# -> 497640

python3 src/hypercatalan.py --table 3 12
# -> prints H(n) for n = 3..12
```

### 2) Build the paper
The LaTeX source is in `paper/main.tex`. It depends only on standard AMS packages:

```bash
cd paper
pdflatex main.tex
# (run twice for references)
```

## Contents

```
.
├── LICENSE
├── README.md
├── paper
│   └── main.tex
├── src
│   └── hypercatalan.py
└── .gitignore
```

## References

- Olarte, I., Santos, F. (2019). *Hypercatalan numbers for the Associahedron and Cyclohedron.* arXiv:1906.05764.
