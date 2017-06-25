.. GlobalBioIm Library documentation master file, created by
   sphinx-quickstart on Sun Jun 25 15:32:10 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GlobalBioIm Library's documentation
===================================

This page contains detailled documentation of each function/class of the Library. The documentation is generated automatically from comments within M-files. It thus constitues the most up-to date documentation of the Library. 

.. default-domain:: mat

Linear Operators (LinOp)
========================

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
