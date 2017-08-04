Abstract classes
****************

This section describes the abstract classes of the GlobalBioIm library. It provides general properties for 
every derived classes.

.. automodule:: Core.Abstract

Map
---

.. autoclass:: Map
    :show-inheritance:
    :members: apply, applyJacobianT, applyInverse, makeComposition

Cost
----

.. autoclass:: Cost
    :show-inheritance:
    :members: apply, applyGrad, applyProx, applyProxFench, makeMapComposition,  applyProx_, applyProxFench_, set_y, set_H
 