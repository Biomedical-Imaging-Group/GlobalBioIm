.. include global.rst

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

