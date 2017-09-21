Important Informations
**********************

This section contains useful informations for developers/users who want to 
   - Implement new classes :class:`Map`, :class:`LinOp`, :class:`Cost` or :class:`Opti`, 
   - Implement new tools in the methods of existing classes (e.g. a fast computation of a gradient of a cost),
   - Use the library for practical problems (see also the provided `examples <http://bigwww.epfl.ch/algorithms/globalbioim/examples.html>`_)

General Philosophy
------------------

The implementation and organization of the library is guided by the following principles:
   - Low code duplication, 
   - Implementing a new map/cost/algorithm should require to edit/create only one file,
   - There should be a solution for implementing new needs (i.e. few hard constraints to respect in implementations),
   - Always have the choice between fast but memory consuming / memory efficient but slower.

Classes Hierarchy
-----------------

The Library is based on 4 abstract classes :class:`Map`,  :class:`LinOp`,  :class:`Cost` and  :class:`Opti` from which 
inherit all the classes used in practice to solve inverse problems. Note that non-linear operators inherit directly from the 
:class:`Map` class since they do not require new specific attributes/methods.

.. figure:: ClassesHierarchy.png
   :scale: 40 %
   :alt: Classes Hierarchy diagram
   :align: center

   Fig 1. Classes Hierarchy diagram of the GlobalBioIm Library.

Interface methods and Core methods
----------------------------------

Classes (except :class:`Opti` classes) contains two types of methods:

   - **Interface methods** which are the methods that can be called from an instanciated object of the class. They perform 
     checks concerning the size conformity of the inputs given to the method and call the corresponding core method. Morever, they
     manage the memoize mechanism. Their implementation is done at the level of the abstract classes. Hence, for an implementation
     of a new class, developpers do not have to deal with the memoize system as well as the checking of inputs sizes.
   - **Core methods** These are the methods that contains implementation. Developpers only need to implement these methods.

Hence, methods always go by pairs (Interface + Core). They have the same name which is followed by an "_" symbol for core ones. 
For instance
   - apply()
   - apply_()


**Note** There is an exception for the :meth:`mtimes` method (i.e`. the overloading of the * operator) which do not have
an associated *core* method. In fact this methods call either :meth:`apply` or :meth:`makeComposition` depending on the 
right hand side of the multiplication.

Memoize and Precomputation options
----------------------------------

The :class:`Map` class provides two attributes which are
    - :attr:`memoizeOpts` a structure of booleans with one field per methods in the class (default all false). If for example 
      the field *memoizeOpts.apply* is set to true then the results of the :meth:`apply` method \\(\\mathrm{y=Hx}\\) is saved.
      Then, if the next call to the :meth:`apply` method is for the same \\(\\mathrm{x}\\), the saved value \\(\\mathrm{y}\\) is directly 
      returned without any computation.
    - :attr:`doPrecomputation` a boolean (default false). When *true*, some methods of the instanciated 
      object will be accelerated at the price of a larger memory consumption. It depends on how the implementation of 
      the class has been done. Hence, if one want to accelerate a method thanks to a precomputation of some quatities, this
      has to be done **when the doPrecomputation is activated** and not by default. This let the possibility to
      avoid the precomputation in case of memory limitation issues.

Let us look at some examples. Consider a convolution oprerator :class:`LinOpConv` 

.. code:: matlab

  H=LinOpConv(fft2(psf));
  H.memoizeOpts.apply=true;  

for a given PSF (\\(512\\times 512 \\times 256 \\)) and for which we activate the memoize option for the :meth:`apply` method.
Then, let make the following calls to :meth:`apply` method.

.. code:: matlab

  >> x=rand(size(psf));
  >> tic;y=H*x;toc;
  Elapsed time is 2.414025 seconds.
  >> tic;y=H*x;toc;
  Elapsed time is 0.085205 seconds.
  >> x(5)=2;
  >> tic;y=H*x;toc;
  Elapsed time is 2.465424 seconds.
  >> tic;y=H*x;toc;
  Elapsed time is 0.083087 seconds.


Compositions of Maps
--------------------

About properties
----------------

Explain about isComplexIn and IsComplexOut

List the "non protected" properties that could actually be modified for really specific reasons but should not in general.

Use the provided templates!
---------------------------

Templates for implementing new :class:`Map`, :class:`LinOp`, :class:`Cost` or :class:`Opti` are provided to help developers.
They can be found under the names:
 - TemplateMap.m
 - TemplateLinOp.m
 - TemplateCost.m
 - TemplateOpti.m