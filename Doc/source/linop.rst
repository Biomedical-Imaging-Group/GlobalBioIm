Linear Operators (LinOp)
************************

This section contains linear operator classes which all derive from the abstract class :class:`LinOp`.

.. warning:: Some methods defined in the main class :class:`LinOp` are **abstract** which means that they need to be implemented in the derived classes. However, all the derived linear operators do not implement all these abstract methods. For example, if one implements a non-invertible linear operator, the method :meth:`inverse` do not have any sense and is thus not implemented for this linear operator. Hence, in the following, abstract methods which are **not mentionned** into the **documentation** of a derived linear operator are **not implemented** for this linear operator. Concerning the non-abstract methods of class :class:`LinOp`, they are always inherited by all subclasses and are not mentionned in their respective documentation (except if there is a reimplementation for any reasons).

.. automodule:: LinOp

LinOpIdentity
-------------

.. autoclass:: LinOpIdentity
    :show-inheritance:
    :members: apply, applyAdjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 
    
LinOpGrad
---------

.. autoclass:: LinOpGrad
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, makeHHt

LinOpHess
---------

.. autoclass:: LinOpHess
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 

LinOpSum
--------

.. autoclass:: LinOpSum
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 

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