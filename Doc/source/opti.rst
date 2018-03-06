.. include global.rst

Optimization Algorithms (Opti)
==============================

	This section contains optimization algorithms classes which all derive from the abstract class :class:`Opti`.
	
.. automodule:: Opti

OptiADMM
--------

.. autoclass:: OptiADMM
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, initialize, doIteration, updateParams

OptiChambPock
-------------

.. autoclass:: OptiChambPock
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, initialize, doIteration, updateParams

OptiConjGrad
------------

.. autoclass:: OptiConjGrad
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, initialize, doIteration, updateParams


OptiDouglasRachford
-------------------

.. autoclass:: OptiDouglasRachford
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, initialize, doIteration, updateParams

OptiFBS
-------

.. autoclass:: OptiFBS
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, updateSet, initialize, doIteration, updateParams

OptiFGP
-------

.. autoclass:: OptiFGP
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, updateSet, initialize, doIteration, updateParams, setBounds, setLambda
    
OptiGradDsct
------------

.. autoclass:: OptiGradDsct
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, set_b, initialize, doIteration, updateParams

OptiPrimalDualCondat
--------------------

.. autoclass:: OptiPrimalDualCondat
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, initialize, doIteration, updateParams

OptiRichLucy
------------

.. autoclass:: OptiRichLucy
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, initialize, doIteration, updateParams

OptiVMLMB
---------

.. autoclass:: OptiVMLMB
    :show-inheritance:
    :members: run, starting_verb, ending_verb, test_convergence, initialize, doIteration, updateParams

OutputOpti
----------

.. autoclass:: OutputOpti
    :show-inheritance:
    :members: init, update

.. automodule:: Opti.TestCvg

TestCvg
-------

.. autoclass:: TestCvg 
    :show-inheritance:   
    :members: init, update
