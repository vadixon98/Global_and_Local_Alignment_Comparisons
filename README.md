# Sequence Alignment

A suite of two R scripts for performing **global** and **local** sequence alignments using dynamic programming, complete with graphical edit-graph visualizations and console-based alignment printing. 

---

## Requirements

* **R** (version 4.0 or higher)
* No external packages required (uses base R functions and `graphics` for plotting).

---

## Installation

1. Download or clone the repository:

   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```
2. Ensure R is installed and accessible via your command line or RStudio.

---

## Scripts Overview

This repository contains the following two R scripts:

### 1. `Find Global Sequence Alignment.R`

Performs Needleman–Wunsch global alignment:

* **Functions**:

  * `find_gsa(v, w, sigma, mu)`: Constructs the global score matrix (`s`) and backtrace pointer matrix (`b`) for sequences `v` and `w` using gap penalty `sigma` and mismatch penalty `mu`.
  * `print_alignment(b, v, w, i, j)`: Recursively follows the backtrace pointers to print aligned sequence pairs (with `"-"` representing gaps) to the console.

### 2. `Global and Local alignment comparisons with graphs.R`

Combines both local (Smith–Waterman) and global alignment routines with graphical output:

* **Predefined Data**:

  * Short example sequences `s1` and `s2`.
  * Long haplotype sequences `s5` and `s6`.
* **Functions**:

  * `find_lsa(v, w, sigma, mu)`: Computes local alignment score matrix (`s`) and backtrace pointers (`b`), resetting to zero for negative scores.
  * `find_gsa(v, w, sigma, mu)`: Same global alignment as in script #1.
  * `plot_edit_graph(dst, backtrace)`: Plots an edit graph, marking match/mismatch cells and highlighting maximum-scoring positions.
  * `print_local_alignment(q)`: Prints optimal local alignment (from `find_lsa()`) to the console.
  * `gather_gsa(b, v, w, i, j, fn)` and `print_gsa(b, v, w, i, j, fn)`: Reconstruct and write global alignments to a file in a two-line FASTA-like format.

---

## Usage Example

```r
# 1. Global alignment only
source("Find Global Sequence Alignment.R")

# Perform global alignment of sequences 's1' and 's2'
g <- find_gsa(s1, s2, sigma = 2, mu = 1)
# Print to console
print_alignment(g$b, g$vectorV, g$vectorW, nrow(g$b), ncol(g$b))

# 2. Combined global & local with graphs
source("Global and Local alignment comparisons with graphs.R")

# Local alignment and graph
q <- find_lsa(s1, s2, sigma = 2, mu = 1)
plot_edit_graph(q$s, q$b)
print_local_alignment(q)

# Global alignment and file export
g2 <- find_gsa(s5, s6, sigma = 2, mu = 1)
print_gsa(g2$b, g2$vectorV, g2$vectorW, nrow(g2$b), ncol(g2$b), fn = "global_alignment_output.txt")
```

---

## Usage Tips

* **Gap & Mismatch Penalties**: Choose `sigma` (gap) and `mu` (mismatch) values according to your substitution scoring needs.
* **Memory**: Both scripts use `O(n*m)` memory; very long sequences may require substantial RAM.
* **Visualization**: Use `plot_edit_graph()` to inspect high-scoring regions before printing alignments.
* **Console vs. File Output**:

  * Use `print_alignment()` and `print_local_alignment()` for quick on-screen checks.
  * Use `print_gsa()` to generate reusable alignment files.

---

## Contributing

Contributions, bug reports, and enhancements are welcome. Please fork the repository, commit your changes, and submit a pull request.

---

## License

This project is released under the MIT License. Include a copy of the license when distributing.
