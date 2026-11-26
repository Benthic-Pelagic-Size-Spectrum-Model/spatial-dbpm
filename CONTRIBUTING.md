# Contributing to spatial-dbpm (dbpmr)

Thanks for contributing to the `dbpmr` benthic–pelagic size-spectrum model. A few
conditions and a simple workflow keep the model coherent and properly credited.

## Conditions of use

- **Work stays in this organisation.** All DBPM development belongs in the
  Benthic-Pelagic-Size-Spectrum-Model organisation. Please contribute changes back
  here via pull request rather than maintaining a diverging public fork.
- **Co-authorship.** Analyses, papers, and other outputs that use the model are
  expected to offer co-authorship to the model developers. Please get in touch
  early about intended use.
- **Licence.** The code is GPL-3. By contributing, you agree your contributions are
  licensed under GPL-3.

## Workflow

`main` is protected: every change lands through a reviewed pull request — no direct
pushes to `main`, for anyone.

1. **Branch** off `main` (collaborators) or **fork** (external contributors). You do
   **not** need elevated permissions — `Write` is enough to push a branch and open a PR.
2. Make a focused change on a feature branch.
3. **Open a pull request** against `main`, describing what changed and why.
4. A **code owner reviews** (see [`.github/CODEOWNERS`](.github/CODEOWNERS)); address
   comments, then it is squash-merged once approved.

## Expectations for engine changes

- The C engine is `dbpmr/src/SizeSpectra.c`. Changes there get careful maintainer review.
- Add or update **tests** (`dbpmr/tests/testthat/`) for any behavioural change.
- **Preserve the 0-D baseline:** a `spatial_dim = 0` run must stay bit-identical unless
  the change is deliberately to that path. New engine toggles default to off / identity.
- Rebuild cleanly before pushing: `rm -f dbpmr/src/*.o dbpmr/src/*.so`, then run the
  test suite (`devtools::test()`).

## Questions

Open an issue, or contact the maintainer (Julia Blanchard).
