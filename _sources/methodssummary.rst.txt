List of Methods
***************

This section gives an overview of all the methods for each abstract class (and thus all inherited classes). All the 
core methods (i.e. with "_") can be reimplemented in any derived class (if relevant). Otherwise the default implementation
of the mother class is used.

**Note** the :meth:`applyJacobianT` method of :class:`Map` becomes sealed (i.e. not re-implementable in derived classes) from
the level of :class:`LinOp` and :class:`Cost`. The reason is that for these classes, the Jacobian is respectively the
adjoint and the gradient. Instead, the methods :meth:`applyAdjoint` and :meth:`applyGrad` are implemented.

Map
---

+-------------------------+--------------------+------------------------------------------------------+
| Interface Methods       | Core Methods       | Description                                          |
+=========================+====================+======================================================+
| apply()                 | apply_()           | | Apply the :class:`Map` to the given input x        |
+-------------------------+--------------------+------------------------------------------------------+
| applyJacobianT()        | applyJacobianT_()  | | Apply the transpose of the Jacobian operator at x  |
|                         |                    | | to a given vector v.                               |
+-------------------------+--------------------+------------------------------------------------------+
| applyInverse()          | applyInverse_()    | | Apply the inverse :class:`Map` to the given input x|                                                      
+-------------------------+--------------------+------------------------------------------------------+
| makeInversion()         | makeInversion_()   | | Build the inverse :class:`Map`                     |  
+-------------------------+--------------------+------------------------------------------------------+
| makeComposition()       | makeComposition_() | | Compose a :class:`Map` with another :class:`Map`.  |                                                   
+-------------------------+--------------------+------------------------------------------------------+
| plus()                  | plus_()            | | overload + operator for Maps.                      |
+-------------------------+--------------------+------------------------------------------------------+
| minus()                 | minus_()           | | overload - operator for Maps.                      |
+-------------------------+--------------------+------------------------------------------------------+
| mpower()                | mpower_()          | | overload ^ operator for Maps (i.e. M.^-1 build the |
|                         |                    | | inverse :class:`Map`).                             |
+-------------------------+--------------------+------------------------------------------------------+
| mtimes()                | N/A                | | overload * operator for Maps. If the input x is    |
|                         |                    | | an array, calls apply(x). If it is a :class:`Map`  |
|                         |                    | | object, calls makeComposition(x).                  |
+-------------------------+--------------------+------------------------------------------------------+
| times()                 | times_()           | | overload .* operator for Maps.                     |
+-------------------------+--------------------+------------------------------------------------------+
| size()                  | N/A                | | overload size function for Maps.                   |
+-------------------------+--------------------+------------------------------------------------------+

LinOp
-----

+-------------------------+-----------------------+------------------------------------------------------+
| Interface Methods       | Core Methods          | Description                                          |
+=========================+=======================+======================================================+
| apply()                 | apply_()              | | Inherited from :class:`Map`                        |
+-------------------------+-----------------------+------------------------------------------------------+
| applyAdjoint()          | applyAdjoint_()       | | Apply the adjoint operator to point x              |
+-------------------------+-----------------------+------------------------------------------------------+
| applyHtH()              | applyHtH_()           | | Apply the HtH operator to point x                  |
+-------------------------+-----------------------+------------------------------------------------------+
| applyHHt()              | applyHHt_()           | | Apply the HHt operator to point x                  |
+-------------------------+-----------------------+------------------------------------------------------+
| applyAdjointInverse()   | applyAdjointInverse_()| | Apply the inverse of adjoint operator to point x   |
+-------------------------+-----------------------+------------------------------------------------------+
| makeAdjoint()           | makeAdjoint_()        | | Build a :class:`LinOp` that implements the adjoint.|
+-------------------------+-----------------------+------------------------------------------------------+
| makeHtH()               | makeHtH_()            | | Build a :class:`LinOp` that implements HtH.        |
+-------------------------+-----------------------+------------------------------------------------------+
| makeHHt()               | makeHHt_()            | | Build a :class:`LinOp` that implements HHt.        |
+-------------------------+-----------------------+------------------------------------------------------+
| applyInverse()          | applyInverse_()       | | Inherited from :class:`Map`                        |                                                      
+-------------------------+-----------------------+------------------------------------------------------+
| makeInversion()         | makeInversion_()      | | Inherited from :class:`Map`                        |  
+-------------------------+-----------------------+------------------------------------------------------+
| makeComposition()       | makeComposition_()    | | Inherited from :class:`Map`                        |                                              
+-------------------------+-----------------------+------------------------------------------------------+
| plus()                  | plus_()               | | Inherited from :class:`Map`                        |
+-------------------------+-----------------------+------------------------------------------------------+
| minus()                 | minus_()              | | Inherited from :class:`Map`                        |
+-------------------------+-----------------------+------------------------------------------------------+
| mpower()                | mpower_()             | | Inherited from :class:`Map`                        |
+-------------------------+-----------------------+------------------------------------------------------+
| mtimes()                | N/A                   | | Inherited from :class:`Map`                        |
+-------------------------+-----------------------+------------------------------------------------------+
| times()                 | times_()              | | overload .* operator for Maps.                     |
+-------------------------+-----------------------+------------------------------------------------------+
| size()                  | N/A                   | | Inherited from :class:`Map`                        |
+-------------------------+-----------------------+------------------------------------------------------+

Cost
----

+-------------------------+--------------------+------------------------------------------------------+
| Interface Methods       | Core Methods       | Description                                          |
+=========================+====================+======================================================+
| apply()                 | apply_()           | | Inherited from :class:`Map`                        |
+-------------------------+--------------------+------------------------------------------------------+
| applyGrad()             | applyGrad_()       | | Apply the gradient of the cost to the given x.     |
+-------------------------+--------------------+------------------------------------------------------+
| applyProx()             | applyProx_()       | | Apply the prox of the cost to the given x.         |
+-------------------------+--------------------+------------------------------------------------------+
| applyProxFench()        | applyProxFench_()  | | Apply the prox of the Fenchel transform of the     |
|                         |                    | | cost to the given x.                               |
+-------------------------+--------------------+------------------------------------------------------+
| applyInverse()          | applyInverse_()    | | Inherited from :class:`Map`                        |                                                      
+-------------------------+--------------------+------------------------------------------------------+
| makeComposition()       | makeComposition_() | | Inherited from :class:`Map`                        |                                               
+-------------------------+--------------------+------------------------------------------------------+
| plus()                  | plus_()            | | Inherited from :class:`Map`                        |    
+-------------------------+--------------------+------------------------------------------------------+
| minus()                 | minus_()           | | Inherited from :class:`Map`                        |    
+-------------------------+--------------------+------------------------------------------------------+
| mpower()                | mpower_()          | | Inherited from :class:`Map`                        |    
+-------------------------+--------------------+------------------------------------------------------+
| mtimes()                | N/A                | | Inherited from :class:`Map`                        |   
+-------------------------+--------------------+------------------------------------------------------+
| times()                 | times_()           | | overload .* operator for Maps.                     | 
+-------------------------+--------------------+------------------------------------------------------+
| size()                  | N/A                | | Inherited from :class:`Map`                        |    
+-------------------------+--------------------+------------------------------------------------------+

Opti
----

+-------------------------+------------------------------------------------------+
| Interface Methods       | Description                                          |
+=========================+======================================================+
| run()                   | | Run the algorithm from a given initial point.      |
+-------------------------+------------------------------------------------------+
| initialize()            | | Initialize the algorithm (e.g. auxilliary var).    |
+-------------------------+------------------------------------------------------+
| doIteration()           | | Performs one iteration of the algorithm.           |
+-------------------------+------------------------------------------------------+
| updateParams()          | | Update algorithm parameters (e.g. descent step).   |
+-------------------------+------------------------------------------------------+
| starting_verb()         | | Display starting message.                          |
+-------------------------+------------------------------------------------------+
| ending_verb()           | | Display ending message.                            |
+-------------------------+------------------------------------------------------+
