# Reporting an issue

Open an [issue](https://github.com/byuzbasi/emaxProjection/issues) with a minimal reproducible example,
or contact bahadir.yuzbasi@inonu.edu.tr. Include the package
version, `sessionInfo()`, operating system, compiler information when relevant,
the dose-axis units, calibration, tolerance, split budget, and complete numerical
status. Use synthetic or properly shareable aggregate data; do not send private
participant records.

Distinguish a package error from an unresolved search or a violation of the
Gaussian common-variance assumptions. Unexpected empty sets must be reported
with their input and status, not silently discarded.

# Proposed changes

Discuss changes to statistical assumptions, kernels, numerical cushions,
defaults, or status conventions with the author before implementation.
Scientific changes require a new validation record and comparison with the
research implementation. Documentation corrections should identify the page
and proposed wording. Discuss proposed changes in the repository issue tracker.
