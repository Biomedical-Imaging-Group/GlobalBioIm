.. GlobalBioIm Library documentation master file, created by
   sphinx-quickstart on Sun Jun 25 15:32:10 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

   
Welcome to the GlobalBioIm Library webpage!
*******************************************

Objectives
----------

When being confronted with a new imaging problem, the common experience is that one has to reimplement 
(if not reinvent) the wheel (=forward model + optimization algorithm), which is very time consuming and 
also acts as a deterrent for engaging in new developments. This Matlab library aims at simplifying this 
process by decomposing the workflow onto smaller modules, including many reusable ones since several aspects
such as regularization and the injection of prior knowledge are rather generic. It also capitalizes on the
strong commonalities between the various image formation models that can be exploited to obtain fast, 
streamlined implementations.

.. figure:: button.png
   :scale: 40 %
   :align: center
   :target: https://c4science.ch/diffusion/2843/


This page contains detailled documentation of each function/class of the Library. The documentation is generated 
automatically from comments within M-files. It thus constitues the most up-to date documentation of the Library. 


News
----

  - June 2017: First public release of the library (v0.1) `get it <https://c4science.ch/diffusion/2843/>`_!

Reference
---------

  - M. Unser, E. Soubies, F. Soulez, M. McCann, L. Donati, 
    `GlobalBioIm: A Unifying Computational Framework for Solving Inverse Problems <http://bigwww.epfl.ch/publications/unser1701.html>`_ 
    Proceedings of the OSA Imaging and Applied Optics Congress on Computational Optical Sensing and Imaging (COSI'17), San Francisco CA, USA, June 26-29, 2017, paper no. CTu1B.

Contents
--------
 
.. toctree::
   :maxdepth: 1
 
   abstract
   linop
   nonlinop
   cost
   opti
   infos
   methodssummary
   examples
   Download/Clone (v.01) <https://c4science.ch/diffusion/2843/>


 
.. default-domain:: mat





