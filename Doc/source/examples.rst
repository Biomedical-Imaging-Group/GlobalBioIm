
.. _ref-examples:

Examples
********

This section presents some examples of use of the Library. 

3D Deconvolution
----------------

We consider the C. elegans embryo real dataset (\\(672 \\times 712 \\times 104 \\)) which can be downloaded `here <http://bigwww.epfl.ch/deconvolution/bio/>`_ and which is presented on Fig. 1. 
This sample contains three kind of structures:
   - chromosomes in the nuclei (DAPI Channel - blue)
   - point-wise spots  (CY3 Channel - red)
   - microtubules (FITC channel - green)
The PSF for each channel, generated from a theoretical model, are also provided `here <http://bigwww.epfl.ch/deconvolution/bio/>`_.
In the following, each channel is deconvolved separately using the same code presented below.

.. figure:: ConvData.png
   :scale: 70 %
   :alt: C. Elegans embryo
   :align: center

   Fig 1. C. Elegans embryo. Three kind of structures: chromosomes in the nuclei (blue), microtubules (green) and a protein 
   stained with CY3 (red).

We start by reading the data 

.. code:: matlab

   %% Reading data
   psf=double(loadtiff(psfname));psf=psf/sum(psf(:));
   y=double(loadtiff(dataname));maxy=max(y(:));y=y/maxy;
   sz=size(y);

We then resize the PSF in order to properly deal with the periodic assumption of the FFT (principally in z where  there is signal at the boundaries of the volume). Note that we use the function *fft_best_dim* (Util/ folder) which allows to find sizes that are
suited for efficient FFT operations.

.. code:: matlab

    padsz=[0,0,52];
    sznew=fft_best_dim(sz+2*padsz);
    halfPad=(sznew-sz)/2;
    psf=padarray(psf,halfPad,0,'both');

We can now define our data fidelity term, TV regularization, and positivity constraint.

.. code:: matlab

    %% Least-Suares Data Fidelity Term
    H=LinOpConv(fftn(fftshift(psf)));                      % Convolution Operator  
    H.memoizeOpts.apply=true;                                         
    S=LinOpSelectorPatch(sznew,halfPad+1,sznew-halfPad);   % Selector Operator
    L2=CostL2(S.sizeout,y);                                % L2 cost function
    LS=L2*S;                                               % Least-Sqaures data term
    %% TV Regularization
    Freg=CostMixNorm21([sznew,3],4);      % TV regularizer: Mixed Norm 2-1
    Opreg=LinOpGrad(sznew);               % TV regularizer: Operator Gradient
    Opreg.memoizeOpts.apply=true;  
    lamb=2e-6;                            % Regularization parameter
    %% Positivity constraint
    pos=CostNonNeg(sznew);                % Non-Negativity: Indicator function
    Id=LinOpIdentity(sznew);              % Identity Operator 

Here, our cost function is as follows
$$ \\mathcal{C}(\\mathrm{x}) = \\frac12 \\|\\mathrm{SHx - y} \\|_2^2 + \\lambda \\|\\nabla \\mathrm{x} \\|_{2,1} + i_{\\geq 0}(\\mathrm{x})$$
where \\(\\lambda >0\\) is the regularization parameter, \\(\\mathrm{S}\\) a selector operator that selects the "non-padded"
part of the convolution result \\(\\mathrm{Hx}\\), and the two others terms are respectively the TV regularization and the positivity
constraint. We are now ready to instanciate and run our algorithm (here ADMM) in order to minimize the above functional.

.. code:: matlab

      maxIt=300;                     % Maximal number of iterations
      Fn={LS,lamb*Freg,pos};         % Functionals F_n constituting the cost 
      Hn={H,Opreg,Id};               % Associated operators H_n
      rho_n=[1e-3,1e-3,1e-3];        % Multipliers rho_n
      ADMM=OptiADMM([],Fn,Hn,rho_n,[],OutputOpti(1,gt,round(maxIt/10)));
      ADMM.ItUpOut=round(maxIt/10);
      ADMM.maxiter=maxIt;
      ADMM.run(xopt);

Here, three splitting have been done: \\(\\mathrm{u_1=Hx}, \\; \\mathrm{u_2=\\nabla x}\\) and \\(\\mathrm{u_3=x}\\).
Note that we do not need to give a solver to the ADMM algorithm (5th argument) since the library operator algebra makes that
building the operator
$$\\rho_1 \\mathrm{H^*H} + \\rho_2 \\nabla^* \\nabla + \\rho_3 \\mathrm{I}$$
results in a :class:`LinOpConv` which is invertible. Hence ADMM builds this operator automatically and uses its inverse for the
linear step  of the algorithm (minimization over \\(\\mathrm{x}\\)).

The deconvolved image is presented in Fig. 2.

.. figure:: Deconv3D.png
   :scale: 70 %
   :alt: Deconvolution result.
   :align: center

   Fig 2. Deconvolution result of the C. Elegans embryo.
