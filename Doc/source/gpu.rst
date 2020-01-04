.. _ref-gpu:

Speedup with GPU
****************

Main functions
--------------

Any script implemented using the library can be easily switched to GPU computation. 
To that end, one only needs two simple functions:
   - useGPU
   - gpuCpuConverter

The function :meth:`useGPU` takes one single argument which can be

   - 0 : computation is done on CPU (default)
   - 1 : computation is done on GPU using the `Matlab Parrallel Computing Toolbox <https://ch.mathworks.com/help/distcomp/>`_ 
   - 2 : computation is done on GPU using CudaMat (`GitHub <https://github.com/RainerHeintzmann/CudaMat>`_, `Documentation <http://www.nanoimaging.de/CudaMat/>`_) 



The function :meth:`gpuCpuConverter` allows to convert numeric variables to the correct type according to the option selected with 
:meth:`useGPU` (respectively, double, gpuarray or cuda). For instance, if x is a :class:`double` and :attr:`useGPU` is set to 0, then :attr:`x=gpuCpuConverter(x)` will not change x. 
But if x is a :class:`double` and :attr:`useGPU` is set to 1, then :attr:`x=gpuCpuConverter(x)` will convert x to a :class:`gpuarray` 
(and send it to the graphic card). Hence, for each created/loaded variable "x" which is sufficiently large
(i.e. not for small vectors), one has to add the following line to the script

.. code:: matlab

    x=gpuCpuConverter(x);

Then, the switch between CPU and GPU computation is simply controlled by :meth:`useGPU` which can be placed at the beginning of the script.


Functions which generate data
-----------------------------

Matlab functions such as :meth:`ones` and :meth:`zeros` require a memory allocation. When GPU is activated, it is faster to generate directly
these data on the graphic card instead of allocating CPU memory and then transfering them to the GPU. In order to make things transparents and
keep the code clear, all the occurences of :meth:`ones` and :meth:`zeros` (i.e. within scripts but also within library classes) should be replaced by 
:meth:`ones_` and :meth:`zeros_`. These two functions have been defined in order to allocate the memory on the CPU or directly on the GPU according to the state of :attr:`useGPU`.


Concrete example
----------------
The 3D deconvolution example provided in the "Example/" folder of the library shows a concrete use of the GPU functionality.
