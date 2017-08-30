.. include global.rst

Cost Functions (Cost)
=====================

	This section contains cost functions classes which all derive from the abstract class :class:`Cost`.
	
.. automodule:: Cost


CostL2
------
.. autoclass:: CostL2
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y
       
CostL1
------
.. autoclass:: CostL1
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostKullLeib
------------
    
.. autoclass:: CostKullLeib
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostMixNorm12
-------------
    
.. autoclass:: CostMixNorm12
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostMixNormSchatt1
------------------
    
.. autoclass:: CostMixNormSchatt1
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y
    
CostComplexHyperBolic
---------------------
    
.. autoclass:: CostComplexHyperBolic
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y