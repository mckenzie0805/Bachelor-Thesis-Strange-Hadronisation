# Yields vs Multiplicity

This directory contains ROOT macros used to study the dependence of strange hadron yields on charged-particle multiplicity in proton–proton collisions. The analysis focuses on relative yield trends and comparisons between different PYTHIA 8 hadronisation models.

## Files

- `avgyieldstopions_vs_multiplicity.C` – Computes and plots average strange hadron yields normalised to pion yields as a function of event multiplicity, enabling direct comparison between Monash and junction-enhanced models.

- `multiplicity_distribution.C` – Produces charged-particle multiplicity distributions used to define multiplicity classes for yield-based analyses.

## Notes

- Yield normalisation to pions is used to reduce trivial multiplicity scaling effects.
- The analysis compares PYTHIA 8 Monash and junction-enhanced colour reconnection models.
- Experimental ALICE data are not directly overlaid in these plots; the focus is on relative model behaviour.
- These macros support the multiplicity-dependent yield results presented in the Bachelor thesis.

