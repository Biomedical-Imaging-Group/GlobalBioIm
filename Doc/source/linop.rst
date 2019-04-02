Linear Operators (LinOp)
************************

This section contains linear operator classes which all derive from the abstract class :class:`LinOp`.

.. automodule:: LinOp
    
LinOpBroadcast
--------------

.. autoclass:: LinOpBroadcast
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpBroadcastMatrix
--------------------

.. autoclass:: LinOpBroadcastMatrix
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpConv
---------

.. autoclass:: LinOpConv
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpDiag
---------

.. autoclass:: LinOpDiag
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_


LinOpGrad
---------

.. autoclass:: LinOpGrad
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_


LinOpHess
---------

.. autoclass:: LinOpHess
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpMatrix
-----------

.. autoclass:: LinOpMatrix
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpSDFT
---------

.. autoclass:: LinOpSDFT
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_


LinOpShape
----------

.. autoclass:: LinOpShape
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpSum
--------

.. autoclass:: LinOpSum
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

      
LinOpXRay
---------

.. autoclass:: LinOpXRay
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_



.. automodule:: LinOp.SelectorLinOps

SelectorLinOps
--------------

LinOpSelector
.............

.. autoclass:: LinOpSelector
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpDownsample
...............

.. autoclass:: LinOpDownsample
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_

LinOpSelectorPatch
..................

.. autoclass:: LinOpSelectorPatch
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_, makeInversion_
