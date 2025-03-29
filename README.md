# mgcv.taps

`mgcv.taps` is an extension of the `mgcv` package that implements **Test for Arbitrary Parametric Structure (TAPS)** smoothers. It provides a unified framework for constructing and testing structured smooth terms that blend user-defined parametric components with adaptive nonparametric components.

## Installation

You can install the development version of `mgcv.taps` from GitHub:

```r
devtools::install_github("harryyiheyang/mgcv.taps")
```

## Overview

This package supports smooth function construction, estimation, and hypothesis testing under arbitrary parametric constraints, by extending the `mgcv` smoothing interface.

Two types of smoothers are currently implemented:

### 1. Univariate Structured Smooth (`bs = "AMatern"`)

This smoother models a function as:

```
f(x) = A(x) %*% alpha + B(x) %*% beta
```

- `A(x)` is a user-defined basis matrix representing a structured parametric component (e.g., `cbind(1, x, x^2)` for quadratic functions, or `cbind(1, x, pmax(0, x - nu))` for piecewise structures).
- `B(x)` is a set of smooth basis functions adaptively constructed using a Matern kernel and constrained to be orthogonal to `A(x)`:

```
t(A) %*% B = 0
```

### 2. Bivariate Structured Smooth (`bs = "A2Matern"`)

This extension handles bivariate smooths with similar structure:

```
f(x1, x2) = A(x1, x2) %*% alpha + B(x1, x2) %*% beta
```

- `A(x1, x2)` can represent structured interaction effects (e.g., `cbind(1, x1, x2, x1 * x2)`).
- `B(x1, x2)` is constructed using Matern kernels and is orthogonal to `A(x1, x2)`.

Additional smoother types may be supported in future versions.

## Estimation

Estimation is performed using `mgcv::gam()`. The `xt` argument should be a list:

- `getA`: a user-specified function returning the parametric basis matrix `A`.
- `para`: optional parameters used in constructing `A`.

### Univariate Example

```r
fit <- gam(y ~ s(x, bs = "AMatern", k = 10,
                 xt = list(getA = function(x, para) cbind(1, x))))
```

### Bivariate Example

```r
fit <- gam(y ~ s(x1, x2, bs = "A2Matern", k = 10,
                 xt = list(getA = function(x1, x2, para) cbind(1, x1, x2, x1 * x2))))
```

## Hypothesis Testing

The package provides a wrapper for Wald tests on specific components of the model using `taps_wald_test`, which actually wraps an unobserved function `mgcv::testStat`:

```r
taps_wald_test(fit, test.component = 1)
```

Alternatively, one can also perform a score test variance component using `taps_wald_test`:

```r
taps_score_test(fit, test.component = 1)
```

Both can be used to formally test the parametric structure encoded in `A`. In these two test functions, `test.component` refers to the index of smooth term to be tested.

## Family of Outcome

Currently, one can use the two novel smoothers, `AMatern` and `A2Matern`, to estimate models under any family of outcome supported by `mgcv`.
For hypothesis testing, `taps_wald_test` also supports all outcome families. However, `taps_score_test` is currently limited to exponential family distributions with canonical link functions.

## Examples

Example usage can be found in the files located in the `example/` directory of the package.

## License

MIT

## Maintainer

Yihe Yang  
Email: yxy1234@case.edu
