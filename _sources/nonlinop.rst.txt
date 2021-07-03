Non-Linear Operators 
********************

This section contains non-linear operator classes which all derive directly from the abstract class :class:`Map`.

.. automodule:: NonLinOp

.. automodule:: NonLinOp.ElementWiseOp

Element-wise Operators
----------------------

OpEWAbs
.......

.. autoclass:: OpEWAbs
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,

OpEWfunc
........

.. autoclass:: OpEWfunc
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,

OpEWInverse
...........

.. autoclass:: OpEWInverse
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,

OpEWSqrt
........

.. autoclass:: OpEWSqrt
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,

OpEWSquaredMagnitude
....................

.. autoclass:: OpEWSquaredMagnitude
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,


