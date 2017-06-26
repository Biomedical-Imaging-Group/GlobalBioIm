.. GlobalBioIm Library documentation master file, created by
   sphinx-quickstart on Sun Jun 25 15:32:10 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GlobalBioIm Library's documentation
===================================

This page contains detailled documentation of each function/class of the Library. The documentation is generated automatically from comments within M-files. It thus constitues the most up-to date documentation of the Library. 

.. note:: This page is under construction, not all the classes available in the library are detailed here yet.

.. default-domain:: mat

Linear Operators (LinOp)
========================

This section contains linear operator classes which all derive from the abstract class :class:`LinOp`.

.. warning:: Some methods defined in the main class :class:`LinOp` are **abstract** which means that they need to be implemented in the derived classes. However, all the derived linear operators do not implement all these abstract methods. For example, if one implements a non-invertible linear operator, the method :meth:`inverse` do not have any sense and is thus not implemented for this linear operator. Hence, in the following, abstract methods which are **not mentionned** into the **documentation** of a derived linear operator are **not implemented** for this linear operator. Concerning the non-abstract methods of class :class:`LinOp`, they are always inherited by all subclasses and are not mentionned in their respective documentation (except if there is a reimplementation for any reasons).

.. automodule:: LinOp

LinOp (abstract class)
----------------------

.. autoclass:: LinOp
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 
    
LinOpIdentity
-------------

.. autoclass:: LinOpIdentity
    :show-inheritance:
    :members: apply, adjoint, HtH, HHt, inverse, adjointinverse, transpose, ctranspose, plus, mtimes, 
    

Cost Functions (Cost)
=====================

	This section contains cost functions classes which all derive from the abstract class :class:`Cost`.
	
.. warning:: Some methods defined in the main class :class:`Cost` are **abstract** which means that they need to be implemented in the derived classes. However, all the derived costs do not implement all these abstract methods. For example, if one implements a non-differentiable cost, the method :meth:`grad` do not have any sense and is thus not implemented for this cost. Hence, in the following, abstract methods which are **not mentionned** into the **documentation** of a derived cost are **not implemented** for this cost. Concerning the non-abstract methods of class :class:`Cost`, they are always inherited by all subclasses and are not mentionned in their respective documentation (except if there is a reimplementation for any reasons).

.. automodule:: Cost

Cost (abstract class)
---------------------

.. autoclass:: Cost
    :show-inheritance:
    :members: eval, grad, eval_grad, prox, prox_fench, o, plus, minus, mtimes
 
CostL2
------
.. autoclass:: CostL2
    :show-inheritance:
    :members: eval, grad, eval_grad, prox, prox_fench, o, plus, minus, mtimes
       
CostL1
------
.. autoclass:: CostL1
    :show-inheritance:
    :members: eval, grad, eval_grad, prox, prox_fench, o, plus, minus, mtimes

CostKullLeib
------------
    
.. autoclass:: CostKullLeib
    :show-inheritance:
    :members: eval, grad, eval_grad, prox, prox_fench, o, plus, minus, mtimes
    
Optimization Algorithms (Opti)
==============================

	This section contains optimization algorithms classes which all derive from the abstract class :class:`Opti`.
	
.. warning:: The :meth:`run` method defined in the main class :class:`Opti` is **abstract** and has to be implemented in any derived class.

.. automodule:: Opti

Opti (abstract class)
---------------------

.. autoclass:: Opti
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence
    
OptiFBS
-------

.. autoclass:: OptiFBS
    :show-inheritance:
    :members: run, 






