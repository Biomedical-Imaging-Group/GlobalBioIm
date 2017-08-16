Important Informations
**********************

This section contains useful informations for developers/users who want to 
   - implement new classes :class:`Map`, :class:`LinOp`, :class:`Cost` or :class:`Opti`, 
   - implement new tools in the methods of existing classes (e.g. a fast computation of a gradient of a cost),
   - use the library for practical problems (see also the provided `examples <http://bigwww.epfl.ch/algorithms/globalbioim/examples.html>`_)

General Philosophy
------------------

Classes Hierarchy
-----------------

Do a diagram with dependances between classes

Interface methods and Core methods
----------------------------------

Explain that only the "underscore" methods have to be implemented (excepted for the special methods, see below). And for use only
the interface methods can be called. Explain the interest of such a structure (check size + memoize)

Special Methods
---------------

Overload of operators +, * ... which do not have Interface + Core methods but only one method
Explain how to overload them in derived classes by calling the superclass method in the default setting

Memoize and Precomputation options
----------------------------------

Explain memoize and precomputation options

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