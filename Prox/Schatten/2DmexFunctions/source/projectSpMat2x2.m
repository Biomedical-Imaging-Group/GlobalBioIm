%Function projectSpMat2x2 finds the following matrix projection
%
% min  ||X-Y||_S(p) 
% ||X||_p <= rho
%
% where Y is a 2x2 real symmetric matrix and p is the order of the corresponding
% Schatten norm.
%
%Matlab Usage: X=projectSpMat2x2(Y,p,rho,c0);
% c0:initial guess for the root.
%
% Function projectSpMat2x2 can handle multiple matrix inputs.
% Y: Y is a multidimensional matrix where its last dimension should be 
% equal to 3 corresponding to the distinct elements of a symmetric 2x2 matrix
% stored column-wise.
% 
% Note that in the multiple matrix input case we can also use a
% different rho for each matrix. Then numel(rho) should be either equal to
% the number of the input matrices or should be a scalar, which means
% that all the matrices are using the same rho.
% =========================================================================
%
%  Author: stamatis.lefkimmiatis@epfl.ch
%
% =========================================================================
