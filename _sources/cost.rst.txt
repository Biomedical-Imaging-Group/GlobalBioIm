.. include global.rst

Cost Functions (Cost)
=====================

	This section contains cost functions classes which all derive from the abstract class :class:`Cost`.
	
.. automodule:: Cost

CostConst
---------
    
.. autoclass:: CostConst
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostGoodRoughness
-----------------
    
.. autoclass:: CostGoodRoughness
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y 

CostHyperBolic
--------------
    
.. autoclass:: CostHyperBolic
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y 

CostKullLeib
------------
    
.. autoclass:: CostKullLeib
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostLinear
----------
.. autoclass:: CostLinear
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostL1
------
.. autoclass:: CostL1
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostL2
------
.. autoclass:: CostL2
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostL2Composition
-----------------
.. autoclass:: CostL2Composition
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostMixNorm21
-------------
    
.. autoclass:: CostMixNorm21
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostMixNorm21NonNeg
-------------------
    
.. autoclass:: CostMixNorm21NonNeg
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostMixNormSchatt1
------------------
    
.. autoclass:: CostMixNormSchatt1
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostRobustPenalization
----------------------
    
.. autoclass:: CostRobustPenalization
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

    
CostTV
------
    
.. autoclass:: CostTV
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y, setProxAlgo

IndicatorFunctions
------------------

.. automodule:: Cost.IndicatorFunctions

CostIndicator
.............
    
.. autoclass:: CostIndicator
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostRectangle
.............

.. autoclass:: CostRectangle
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostReals
.........

.. autoclass:: CostReals
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostNonNeg
..........

.. autoclass:: CostNonNeg
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostComplexRing
...............

.. autoclass:: CostComplexRing
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostComplexCircle
.................

.. autoclass:: CostComplexCircle
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostComplexDisk
...............

.. autoclass:: CostComplexDisk
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y
