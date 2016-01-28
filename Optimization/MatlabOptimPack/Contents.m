% MatlabOptimPack Package 
%
% A Matlab interface to the OptimPack package for constrained or unconstrained 
% optimization using the gradient
%
% The required OptimPack package is available on Eric Thiebaut's home page :
% http://cral.univ-lyon1.fr/labo/perso/eric.thiebaut/?Software/OptimPack
% 
% - m_vmlmb_first.c : MexFile interface to the OptimPack op_vmlmb_first function
% - m_vmlmb_first.m : Help for the Matlab m_vmlmb_first function
%
% - m_vmlmb_next.c : MexFile interface to the OptimPack op_vmlmb_next function
% - m_vmlmb_next.m : Help for the Matlab m_vmlmb_next function
%   WARNING : 
%       To be coherent with the OptimPack library in terms of limited memory 
%       used, but contrarily to the usual Matlab functions, the isave and 
%       dsave variables are passed by reference in the m_vmlmb_next function, 
%       so they might both be considered input and output parameters 
%
% - optim_vmlmb.m : Matlab optimization function using the m_vmlmb_first and 
%                   m_vmlmb_next functions
%
% - test_optim_vmlmb.m : Matlab script to test the package on an unconstrained 
%                        optimization case
% - test_optim_vmlmb.c : C program to test the package on the same unconstrained 
%                        optimization case
%
% - test_optim_vmlmb_const.m  : Matlab script to test the package on an 
%                               unconstrained optimization case
% - test_optim_vmlmb_const.c : C program to test the package on the same 
%                              unconstrained optimization case
%
% - ros.m : Compute the Rosenbrock function used in the Matlab tests scripts
% - grad_ros.m : Compute the gradient of the Rosenbrock function used in the 
%                Matlab tests scripts
%
% Installation :
%	- Install the OptimPack library with the -fPIC option
%		replace : CFLAGS = -O2 -Wall
%		with : CFLAGS = -O2 -Wall -fPIC
%		in the OptimPack Makefile
%	- Compile the Matlab MexFiles with Matlab command mex using the OptimPack 
%         library
%	  	>> mex m_vmlmb_first.c -loptimpack
%	  	>> mex m_vmlmb_next.c -loptimpack
%	- Compile the C tests functions using the OptimPack library
%               gcc -o test_optim_vmlmb.o -c test_optim_vmlmb.c -Wall
%               gcc -o test_optim_vmlmb test_optim_vmlmb.o -lm -loptimpack
%               gcc -o test_optim_vmlmb_const.o -c test_optim_vmlmb_const.c -Wall
%               gcc -o test_optim_vmlmb_const test_optim_vmlmb_const.o -lm -loptimpack
%       - the Matlab scripts and C tests programs should give the same results :
%		for an unconstrained optimization case : 
%			test_optim_vmlmb.m and test_optim_vmlmb.c
%		for a constrained optimization case : 
%			test_optim_vmlmb_const.m and test_optim_vmlmb_const.c
%
%  September 2013
%  IRAP, Observatoire Midi-Pyrénées, Toulouse
%  Hervé Carfantan
