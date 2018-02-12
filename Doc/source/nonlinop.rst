Non-Linear Operators 
********************

This section contains non-linear operator classes which all derive directly from the abstract class :class:`Map`.

.. automodule:: NonLinop

.. automodule:: NonLinop.ElementWiseOp

Element-wise Operators
----------------------

OpEWSquareRoot
..............

.. autoclass:: OpEWSquareRoot
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,

.. autoclass:: OpEWSquaredMagnitude
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,

.. autoclass:: OpEWInverse
    :show-inheritance:
    :members: apply_, applyJacobianT_, applyInverse_, plus_, minus_, mpower_, makeComposition_,