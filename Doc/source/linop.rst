Linear Operators (LinOp)
************************

This section contains linear operator classes which all derive from the abstract class :class:`LinOp`.

.. automodule:: LinOp

LinOpDiag
---------

.. autoclass:: LinOpDiag
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_
    
LinOpGrad
---------

.. autoclass:: LinOpGrad
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_

LinOpHess
---------

.. autoclass:: LinOpHess
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_

LinOpConv
---------

.. autoclass:: LinOpConv
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,
      applyAdjoint_, applyHtH_, applyHHt_, applyAdjointInverse_, makeAdjoint_, makeHtH_, makeHHt_

.. automodule:: LinOp.SelectorLinOps

SelectorLinOps
--------------

LinOpSelector
.............

.. autoclass:: LinOpSelector
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 

LinOpDownsample
...............

.. autoclass:: LinOpDownsample
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 

LinOpSelectorPatch
..................

.. autoclass:: LinOpSelectorPatch
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 