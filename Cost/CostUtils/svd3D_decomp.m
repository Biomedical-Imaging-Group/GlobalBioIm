% function [E,V]=svd3D_decomp(X)
%
% Let X be a NxMxKx6 matrix such that:
%   
%   P_mn = [X(n,m,k,1) X(n,m,k,2) X(n,m,k,3)
%           X(n,m,k,2) X(n,m,k,4) X(n,m,k,5)
%           X(n,m,k,3) X(n,m,k,5) X(n,m,k,6)] 
%           
% is a symmetric matrix. Then the present function computes the eigenvalues
% E(n,m,k,1) E(n,m,k,2) E(n,m,k,3) and the eigenvector 
%           V1 = [V(n,m,k,1) V(n,m,k,2)  V(n,m,k,3)] 
%           V2 = [V(n,m,k,4) V(n,m,k,5)  V(n,m,k,6)] 
%           V2 = [V(n,m,k,5) V(n,m,k,8)  V(n,m,k,9)]  
% Hence the function outputs two matrices E of size NxMxKx3 and V of size NxMxKx9.
%  
%  Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.%
%%
