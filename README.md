# GlobalBioIm 
[![DOI](https://zenodo.org/badge/124297869.svg)](https://zenodo.org/badge/latestdoi/124297869)

This Matlab toolbox provides a set of tools (operators, cost functions, optimization algorithms) for solving inverse problems.

## Getting started 

Download or clone this repository and run the script *setMatlabPath.m* which will update your Matlab path with the library.

## Documentation

A detailled documentation can be found <a href="https://biomedical-imaging-group.github.io/GlobalBioIm/" target="_blank">here</a>.


## v1.0 release notes

#### Major Changes and improvements
- Remove deprecated `LinOpScaling`,
- Improve tests with `CheckLinOp`, `CheckMap` and `TestsSummationLinOp' functions
- implement methods `makeAdjoint` and `makeInversion` in `LinOp`
- more simplifications in `Summation`, `Composition`, `Ìnversion` for `Map` and hence `LinOp`and `Cost`
- New `LinOpBroadcast` linear operator
- Add `ElementWiseOp` non linear operators: `OpEWSquareRoot`, `OpEWInverse` and `OpEWSSquaredMagnitude`.  
- New `CostGoodRoughness` cost function
- Introduce positivity constraint in `CostL1`
- New  `CostTV` cost function
- Rework `Opti` optimization class
- Introduce several `TestCvg` classes to assess convergence of optimization according to different criteria
- Automatic download and compilation of OptimPackLegacy for `OptiVMLMB`
- `OptiConjGrad`no longer takes a `W` parameter
- New ÒptiFGP` (Fast gradient proximal) optimization function



#### Bug fixes and minor changes
- Improve doc
- Implement function `cmpSize` to compare size of vectors
- Faster `CostHyperBolic` cost function- Automatic compilation of `C` files needed for `CostMixNormSchatt1` cost function
- Better examples
- Add functions `GenerateData` and `GenerateData` to generate  star phantom- More options in `LinOpConv'
- Implement `make...` methods in LinOpDFT
- Reworked `OutputOpti` printing class



Contributors
-----------------
Mike Mc Cann, Thomas Debarre, Thanh-an Pham, Emmanuel Soubies and Ferréol Soulez

## Reference

[GlobalBioIm: A Unifying Computational Framework for Solving Inverse Problems](http://bigwww.epfl.ch/publications/unser1701.pdf)  <br />
Proceedings of the OSA Imaging and Applied Optics Congress on Computational Optical Sensing and Imaging (COSI‘17), San Francisco CA, USA, June 26-29, 2017, paper no. CTu1B. <br />
M. Unser, E. Soubies, F. Soulez, M. McCann, L. Donati.
