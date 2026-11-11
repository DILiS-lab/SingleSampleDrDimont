## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The words flagged as possibly misspelled in DESCRIPTION are correct: author names (Hiort,
  Kuijjer), "et al." and the standard terms mRNA and omics.

## Test environments

* local macOS 26.7.1, R 4.5.0, aarch64-apple-darwin20: 0 errors, 0 warnings, 4 notes.
  Besides "New submission", the notes come from the local machine: "unable to verify
  current time", "README.md / NEWS.md cannot be checked without pandoc", and "HTML Tidy
  not recent enough".
* Linux (Docker rocker/r-ver, Ubuntu 24.04, aarch64), R 4.6.1, `--no-manual`:
  0 errors, 0 warnings, 1 note ("New submission").
* win-builder (Windows Server 2022, x86_64), R 4.6.2 and R-devel (2026-11-02 r90827):
  0 errors, 0 warnings, 1 note ("New submission" and the possibly misspelled words above).

## Downstream dependencies

There are currently no downstream dependencies for this package.
