.. include global.rst

Optimization Algorithms (Opti)
==============================

	This section contains optimization algorithms classes which all derive from the abstract class :class:`Opti`.
	
.. automodule:: Opti
    
OptiGradDsct
------------

.. autoclass:: OptiGradDsct
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, set_b(

OptiConjGrad
------------

.. autoclass:: OptiConjGrad
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence

OptiFBS
-------

.. autoclass:: OptiFBS
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence

OptiChambPock
-------------

.. autoclass:: OptiChambPock
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence

OptiPrimalDualCondat
--------------------

.. autoclass:: OptiPrimalDualCondat
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence

OptiADMM
--------

.. autoclass:: OptiADMM
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence

OptiDouglasRachford
-------------------

.. autoclass:: OptiDouglasRachford
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence


OptiRichLucy
------------

.. autoclass:: OptiRichLucy
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence

OptiVMLMB
---------

.. autoclass:: OptiVMLMB
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence

OutputOpti
----------

.. autoclass:: OutputOpti
    :show-inheritance:
    :members: init, update