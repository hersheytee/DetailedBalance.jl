# DetailedBalance

[![Build Status](https://github.com/hersheytee/DetailedBalance.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hersheytee/DetailedBalance.jl/actions/workflows/CI.yml?query=branch%3Amain)

The DetailedBalance.jl package aims to provide a set of functions for fast and easy solar cell detailed balance calculations.

Currently the package can calculate the IV curve for a single-junction solar cell given a specified spectrum file, solar cell temperature, and bandgap or set of bandgaps. See example.jl in src for more.

This package is a work in progress, some of the things in the (very short-termm) pipeline include:
* an expansive suite of unit tests
* extensive documentation
* performance improvements
* extension to multi-junction solar cells
* more examples
* open access to more functions
* JuliaCall wrapper to enable use of the package in Python
* making the package open-source

If stuff is broken, please email me at htalathi7@gmail.com.

