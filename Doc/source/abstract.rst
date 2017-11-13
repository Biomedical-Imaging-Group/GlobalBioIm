Abstract classes
****************

This section describes the abstract classes of the GlobalBioIm library. It provides general properties for 
every derived classes.

.. automodule:: Abstract

Map
---

.. autoclass:: Map
    :show-inheritance:
    :members: apply, applyJacobianT, applyInverse, makeComposition, plus, minus, mpower,
      mtimes, size, apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_

LinOp
-----

.. autoclass:: LinOp
    :show-inheritance:
    :members: applyAdjoint, transpose, ctranspose, applyHtH, applyHHt, applyAdjointInverse, makeHtH, makeHHt,
      apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_

Cost
----

.. autoclass:: Cost
    :show-inheritance:
    :members: applyGrad, applyProx, applyProxFench,   
      apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

Opti
----

.. autoclass:: Opti
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence


.. automodule:: Abstract.Compositions

.. _ref-op-on-Maps:

OperationsOnMaps
----------------

The following classes implement basic operations between Map (LinOp and Cost). 
They are not abstract but generally they do not need to be instanciated. 
They are mainly used inside the methods of the abstract classes :class:`Map`, :class:`LinOp` and
:class:`Cost` for the operator algebra machinery.

.. automodule:: Abstract.OperationsOnMaps


MapComposition
..............

.. autoclass:: MapComposition
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_

MapInversion
............

.. autoclass:: MapInversion
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_

MapSummation
..............

.. autoclass:: MapSummation
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_

LinOpAdjoint
............

.. autoclass:: LinOpAdjoint
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_


LinOpComposition
................

.. autoclass:: LinOpComposition
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_


LinOpInversion
..............

.. autoclass:: LinOpInversion
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_

LinOpSummation
..............

.. autoclass:: LinOpSummation
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_

CostComposition
...............

.. autoclass:: CostComposition
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostMultiplication
..................

.. autoclass:: CostMultiplication
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y

CostSummation
.............

.. autoclass:: CostSummation
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyGrad_, applyProx_, applyProxFench_, set_y