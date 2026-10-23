# Independent reproduction example

These original R projection functions complement the RcppArmadillo package.
From this directory, with dependencies already installed, create a local library
and install the source archive downloaded from the package website:

```sh
mkdir -p local-library
R CMD INSTALL --library=local-library ../emaxProjection_0.1.0.tar.gz
Rscript --vanilla reproduce_clinical.R local-library
```

Place the downloaded archive in the repository root for the command above.
The deterministic script checks the reported clinical interval and agreement
between the original R implementation and the C++ port. Aggregate summaries
represent all 369 released observations from DoseFinding::IBScovars (1.4-1).
Dose codes 0–4 are treated as equally spaced; physical dose units are unavailable.
No participant records, random draws or production simulation are included.

The completed 465,000-dataset study used the R research implementation. This
folder supplies its projection functions and completed summary tables, not the
entire study archive. Additional reproduction material is available from the
author on reasonable request. See the package reproduction guide and LICENSE.
