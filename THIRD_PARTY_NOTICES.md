# Third-party components

`simulation/Greedy/R` and `simulation/Greedy/src` contain only the NSBM
inference components used from noisysbmGGM 0.1.2.3, by Valentin Kilian and
Fanny Villers. Copyright/license remain with their authors; GPL-2 text is
in `simulation/Greedy/LICENSE`. Source algorithms are unchanged. The adapter
provides a private compilation environment and uses the original NIG search
and merge. GGM preprocessing, examples, datasets, vignettes and compiled
artifacts are not redistributed here.

VEM is supplied by the external R package noisySBM (paper version 0.1.4),
not copied into VB. SC is the existing local Sun-Cai-style known-null
implementation; it is not claimed to be a byte-identical upstream release.
The research VB code is kept separate from these comparators. A publication
license for the authors' original research code has not been assigned by
this directory-organization step.
