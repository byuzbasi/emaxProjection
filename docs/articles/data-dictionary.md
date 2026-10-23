# Summary-data dictionary and provenance

## Source and sample

The clinical illustration uses `IBScovars` from **DoseFinding 1.4-1**.
The source is the [CRAN package
archive](https://cran.r-project.org/package=DoseFinding). The analysis
retained all **369 released observations**, with dose counts 71, 78, 75,
72, and 73 at codes 0, 1, 2, 3, and 4. No observation was excluded,
imputed, standardized, or transformed for this analysis. The released
outcome is a baseline-adjusted abdominal-pain response score; larger
values indicate greater improvement according to the source
documentation.

The release contains completers, so this is not an intention-to-treat
reconstruction. Active physical doses and spacing are not supplied. The
equal-spacing code convention is explicit and affects the integrated
target.

## Bundled aggregate files

`clinical_dose_summary.tsv` contains dose-level summaries derived from
the unchanged released sample. `clinical_intervals.tsv` contains the
saved completed interval results. The package contains no
individual-level clinical data.

| File | Field | Meaning |
|:---|:---|:---|
| Dose summary | dose | Ordered dose code, 0 through 4; not a physical dose |
| Dose summary | n | Number of released observations at that dose |
| Dose summary | mean | Sample mean of the released response |
| Dose summary | sd | Sample standard deviation with denominator n - 1 |
| Interval results | method / label | Recorded method identifier and display label |
| Interval results | status | Recorded numerical status |
| Interval results | estimate / qhat | Common fitted target estimate and shape estimate |
| Interval results | lower / upper / width | Reported interval endpoints and width in response units |
| Interval results | calibration | Description of the method’s calibration |

## Pooled variance

Let $`s_j`$ be the dose-specific sample standard deviation and $`n_j`$
its count. The pooled within-dose variance for $`m`$ dose levels is
computed as

``` math
S^2=\frac{\sum_{j=1}^m(n_j-1)s_j^2}{\sum_{j=1}^m n_j-m}.
```

This expression pools pure error under the common-variance assumption;
its denominator is the pure-error degrees of freedom. It reproduces the
variance used by the confidence region, up to the precision of the saved
summary text. It is not estimated from deviations about a fitted Emax
curve.

## Reuse and publication status

The source package and aggregate summaries are available from the
package website and GitHub repository. Source-dataset documentation and
applicable reuse conditions should accompany any future deposition.
