## v1.2 release notes
The v1.2 release contains a new graphical user interface (GUI) plus incremental improvements and bug fixes.

#### Graphical User Interface
The GUI allows to intuitively build tailored reconstruction algorithms with minimal effort. The proficiency of the GlobalBioIm GUI is, once the objective function is defined, to automatically propose algorithms that are adapted to the built cost function. Finally, the MATLAB script corresponding to the pipeline is automatically generated. See more at
https://biomedical-imaging-group.github.io/GlobalBioIm/gui.html

#### Improvements
- improvement to save memory in ADMM
- in Opti classes, previous iterate xold is no longer stored  when it is not necessary
- move (deprecated) SNR stuffs from OutpuOpti to OutputSNR
- add Element-wise function operator: OpEWfunc
- add Element-wise Modulus operator: OpEWAbs
- add CostRobustPenalization (Andrews, Beaton-Tukey, Cauchy, Fair, Huber, Logistic, Talwar-Hinich, Welsch-Dennis)
- add the option to put weights in CostHyperBolic
- add constant cost : CostConst
- add option for Nesterov accelerated gradient descent in OptiGradDesc
- add L infinity cost function
- allow deep copy of objects
- make protected most Map/Linop/Cost properties
- add backtracking rule and reduced-step-size rule in OptiFBS
- add LinOpSelectorPlanes
- add memoize options in CostGoodRoughness (like in CostHyperBolic)
- add attribute first index in LinOpSelector

#### Bug Fixes
- fix bug in the adjoint computation of StackMap
- fix LinOpSum norm computation with multiple index
- fix bug OptiFGP when running on GPU
- fix issue in CostLinear with complex entries
- in OptiConjGrad set_b method : also set CostLinear b to new value
- fix pb to restart algorithms from the last iterate of the previous run

## v1.1.2 release notes
The v1.1.2 release does not contain any major changes but mainly incremental improvements and bug fixes. 

#### Improvements 
- New `LinOpXRay` for 2D radon transform (Note that the adjoint is not really the adjoint, but only an approximation (up to about 20 dB),
- LinOp compatibility with sparse matrices,
- `Util/estimateNorm.m` for estimation of the norm of an operator,
- New `OutputOptiSNR` to print SNR during optimization,
- Better summation and multiplication of Map objects with scalars,
- Rename `OpEWSquareRoot` to` OpEWSqrt`,
- Add some kind of unitary tests,
- `CostHyperBolic`: Epsilon can be a vector,
- non scalar step size in `OptiFBS`


#### Bug fixes
Many bug fixes, among them:
- Generalize Compiler option for C/C++ files and reactivate openmp for the Schatten norm
- Better handling of singletons dimensions.



## v1.1.1 release notes
The v1.1.1 release does not contain any major changes but mainly incremental improvements and bug fixes. 

#### Improvements 
- Sfft functions: faster implementation, work with any dimension, new padding option argument.
- Merge LinOpDFT into LinOpSDFT. LinOpDFT remains now as an alias.
- Remove code duplication in SelectorLinOps.
- Add unitary tests (using the provided examples, see /Util/UnitTest/Script.m)
- Add subgradient in CostL1 
- Add test for applyAdjointInverse method in checkMap
- Add property for real input in CostComposition. When activated, takes the real part after gradient computation.

#### Bug fixes
- Fix bug in applyHtH and applyHHt of LinOpSummation
- Fix a memory issue in Cost
- Fix bug in the computation of the Lipschitz constant of CostSummation.
- Fix index of the coefficients vetor in the plus method of LinOp
- Fix superclass constructor call in CostReals

## v1.1 release notes
The v1.1 release contains a new functionality for GPU computation plus incremental improvements and bug fixes. 

#### GPU functionality
- add the functions ones_, zeros_, useGPU, isGPU, gpuCpuConverter
- replace occurences of all() by sum( test)==0 
- replace occurences of ones() and zeros() by  ones_() and zeros_()
- update the 3D deconvolution example with GPU
- add a page in the documentation dedicated to GPU

#### Other changes
- fix a bug in OptiVMLMB compilation (you should recompile it by running installOptimPack.m)
- set saveXopt=false by default in OutputOpti
- add flag isSeparable for Cost classes. Generic implementation of applyProx_ for the sum of a separable cost with an indicator cost.
- add a direct prox computation for the composition of a CostL2 with a Downsampled Convolution
- add CostMixNorm21 a non negativity constrained version of CostMixNorm21
- ADMM can now take a Sum of CostL2Composition in the F0 argument
- remove OutputOptiADMM


## v1.0.1 release notes
The v1.0.1 release does not contain any major changes but mainly incremental improvements and bug fixes. 

#### LinOp
- New LinOpBroadcastMatrix
- fix pb in makeComposition of LinOpSum other cases in makeComposition btw LinOpDiag and LinOpConv
- LinOpHess is now working for any dimensions, with the choice of the index along which the hessian is computed (Like in LinOpGrad). The method makeHtH return the correct corresponding convolution when the boundary condition is circular
- Adapt MixNormSchatten to the new Hessian functionalities
- makeHtH in LinOpGrad works now for any number of dimensions. It also returns the correct LinOpConv if the gradient is applied to specific dimensions
- LinOpCpx is deprecated as it is not really a LinOp

#### Cost
- New CostLinear
- CostL2 now accept any linop as weighting (metric) operator
- Indicator functions: Move the size parameter at the first place according to what was done for all maps

#### Opti
- add methods computeCost and computeSNR in OutputOpti to facillitate deriving classes
- Fix copy on write side effect in OptiVMLMB
- Fix OptiVMLMB such that the returned value is the last accepted one
- Add OutputOptiADMM which evaluates the Lagrangian instead of the cost function
- OutputOpti can deal with both logical or scalar flags

#### Bug fixes other changes
- Improve doc
- Improve examples
- Faster operation when memoize is off
- Add Conjugate Gradient in Deconv_LS_NoReg


## v1.0 release notes

#### Major Changes and improvements
- Remove deprecated `LinOpScaling`,
- Improve tests with `CheckLinOp`, `CheckMap` and `TestsSummationLinOp' functions
- implement methods `makeAdjoint` and `makeInversion` in `LinOp`
- more simplifications in `Summation`, `Composition`, `??nversion` for `Map` and hence `LinOp`and `Cost`
- New `LinOpBroadcast` linear operator
- Add `ElementWiseOp` non linear operators: `OpEWSquareRoot`, `OpEWInverse` and `OpEWSSquaredMagnitude`.  
- New `CostGoodRoughness` cost function
- Introduce positivity constraint in `CostL1`
- New  `CostTV` cost function
- Rework `Opti` optimization class
- Introduce several `TestCvg` classes to assess convergence of optimization according to different criteria
- Automatic download and compilation of OptimPackLegacy for `OptiVMLMB`
- `OptiConjGrad`no longer takes a `W` parameter
- New `OptiFGP` (Fast gradient proximal) optimization function



#### Bug fixes and minor changes
- Improve doc
- Implement function `cmpSize` to compare size of vectors
- Faster `CostHyperBolic` cost function
- Automatic compilation of `C` files needed for `CostMixNormSchatt1` cost function
- Better examples
- Add functions `GenerateData` and `GenerateData` to generate  star phantom- More options in `LinOpConv'
- Implement `make...` methods in LinOpDFT
- Reworked `OutputOpti` printing class


