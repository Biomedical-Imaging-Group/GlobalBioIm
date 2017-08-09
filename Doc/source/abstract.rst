Abstract classes
****************

This section describes the abstract classes of the GlobalBioIm library. It provides general properties for 
every derived classes.

.. automodule:: Abstract

Map
---

.. autoclass:: Map
    :show-inheritance:
    :members: apply, applyJacobianT, applyInverse, makeComposition, plus, minus, 
      mtimes, apply_, applyJacobianT_, applyInverse_, makeComposition_, plus_, minus_

LinOp
-----

.. autoclass:: LinOp
    :show-inheritance:
    :members: applyAdjoint, transpose, ctranspose, applyHtH, applyHHt

Cost
----

.. autoclass:: Cost
    :show-inheritance:
    :members: applyGrad, applyProx, applyProxFench,   
      applyGrad_, applyProx_, applyProxFench_, makeComposition_, set_y


.. automodule:: Abstract.Compositions

OperationsOnMaps
----------------

    The following classes implement basic operations between Map (LinOp and Cost).

.. automodule:: Abstract.OperationsOnMaps

MapComposition
..............

.. autoclass:: MapSummation
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, makeComposition_, plus_, minus_

.. autoclass:: MapComposition
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, makeComposition_, plus_, minus_
 