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

.. automodule:: Opti.OutputOpti

OutputOpti
----------

OutputOpti (Default)
....................

.. autoclass:: OutputOpti
    :show-inheritance:
    :members: init, update, computeCost, computeSNR

OutputOptiConjGrad
..................

.. autoclass:: OutputOptiConjGrad
    :show-inheritance:
    :members: init, update, computeCost, computeSNR

.. automodule:: Opti.TestCvg

TestCvg
-------

TestCvg (Default)
..................

.. autoclass:: TestCvg 
    :show-inheritance:   
    :members: testConvergence

TestCvgCombine
..............

.. autoclass:: TestCvgCombine 
    :show-inheritance:   
    :members: testConvergence

TestCvgCostAbsolute
...................

.. autoclass:: TestCvgCostAbsolute 
    :show-inheritance:   
    :members: testConvergence

TestCvgCostRelative
...................

.. autoclass:: TestCvgCostRelative 
    :show-inheritance:   
    :members: testConvergence

TestCvgStepRelative
...................

.. autoclass:: TestCvgStepRelative 
    :show-inheritance:   
    :members: testConvergence

TestCvgMaxSnr
.............

.. autoclass:: TestCvgMaxSnr 
    :show-inheritance:   
    :members: testConvergence

TestCvgADMM
...........

.. autoclass:: TestCvgADMM 
    :show-inheritance:   
    :members: testConvergence
