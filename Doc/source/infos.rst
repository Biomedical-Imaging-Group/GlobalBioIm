Important Informations
**********************

This section contains useful informations for developers/users who want to 
   - Implement new classes :class:`Map`, :class:`LinOp`, :class:`Cost` or :class:`Opti`, 
   - Implement new tools in the methods of existing classes (e.g. a fast computation of a gradient of a cost),
   - Use the library for practical problems (see also the provided :ref:`examples <ref-examples>`)

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
Then, let us make the following calls to :meth:`apply` method.

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

Here, one can appreciate the effect of this *memoize* option. To observe the effect of the *precomputation* option, 
we instantiate a :class:`CostL2` object which we combine with our convolution operator and for which we activate the 
*precomputation* option.

.. code:: matlab

  >> y=rand(size(psf));  
  >> LS=CostL2([],y); 
  >> F=LS*H;
  >> F.doPrecomputation=1;

Let us evaluate the gradient of the cost *F* at a point x.

.. code:: matlab

  >> x=rand(size(psf));
  >> tic;g=F.applyGrad(x);toc;
  Elapsed time is 5.236043 seconds.
  >> tic;g=F.applyGrad(x);toc;
  Elapsed time is 2.554012 seconds.
  >> tic;g=F.applyGrad(x);toc;
  Elapsed time is 2.572284 seconds.

For a :class:`CostL2`, when the *precomputation* option is activated, the gradient is computed using

$$ \\nabla F(\\mathrm{x})= \\mathrm{H^*Hx - H^*y},$$

which allows to take benefit from a fast implementation of \\( \\mathrm{H^*H}\\) (for the above example \\( \\mathrm{H^*H}\\) is also a 
convolution). Here, at the first call of :meth:`applyGrad`, the quantity \\(\\mathrm{H^*y}\\) is computed and stored 
(hence 4 FFT/IFFT are performed). Then, for all the following calls to :meth:`applyGrad` the compuation is now reduced to
the application of \\( \\mathrm{H^*H}\\) which requires only 2 FFT/IFFT in this case.

Note that in the above example we computed the gradient 3 times over the same x without activating the *memoize* option 
for the :meth:`applyGrad` method in order to show the effect of *precomputation*. Of course, doing new calls to :meth:`applyGrad`
with the same x after having activated the *memoize* option produces,

.. code:: matlab

  >> F.memoizeOpts.applyGrad=true;
  >> tic;g=F.applyGrad(x);toc;
  Elapsed time is 2.572987 seconds.
  >> tic;g=F.applyGrad(x);toc;
  Elapsed time is 0.075061 seconds.



Compositions
------------

The library enjoys a nice operator algebra mechanism that allows sone generic implementations. This is made possible
thank to the methods prefixed by *make* (i.e. :meth:`makeComposition_`, :meth:`makeAdjoint_`, :meth:`makeHtH_`, :meth:`makeHHt_`...).
as well as the :meth:`plus_`, :meth:`minus_` and :meth:`mpower_` methods. By default these methods will instanciate 
some :ref:`Operations on Maps <ref-op-on-Maps>` objects which may lose some properties such as


Auxiliary Utilities
-------------------

The library contains a folder *Util/* with several function. These include viewvers or checking functions which give a
first control that an implementation is correct. For instance, the **CheckMap** function verifies some basic relations
between the different methods implemented in the given :class:`Map`.

.. code:: matlab

    >> H=LinOpConv(fft2(psf));
    >> checkMap(H)
    -- Checking Map with name LinOpConv--
    apply OK
    applyJacobianT OK
    applyInverse OK
        accurate to 2.952430e+02 dB, OK
    -- LinOp-specific checks --
    applyAdjoint OK
        accurate to 2.735641e+02 dB, OK
    applyHtH OK
        accurate to 3.304399e+02 dB, OK
    applyHHt OK
        accurate to 3.167828e+02 dB, OK


Use the provided templates!
---------------------------

Templates for implementing new :class:`Map`, :class:`LinOp`, :class:`Cost` or :class:`Opti` are provided to help developers.
They can be found under the names:
 - TemplateMap.m
 - TemplateLinOp.m
 - TemplateCost.m
 - TemplateOpti.m