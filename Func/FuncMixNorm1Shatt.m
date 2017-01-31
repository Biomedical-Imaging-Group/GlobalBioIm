classdef FuncMixNorm1Shatt < Func
    %% FuncMixNorm1Shatt : Mixed l1-Shatten Norm
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % This function implements the mixed l1-Shatten norm as proposed in [1]:
    %     sum_n || (Hx)(n,:) ||_{Sp},   p=>1 (p=Inf is possible)
    % where H is a LINOP, n can stand for multiple indexes (and ':' for the remaining ones) and Sp
    % denotes the p-order Shatten norm defined by 
    %      ||X||_{Sp} = [sum_k  (sig_k(X))^p]^{1/p},
    % where sig_k(X) is the k-th singular value of X. In other words it is the lp-norm 
    % of the signular values of X.
    %
    % Note: The actual implementation ...
    %
    % -- Example
    % F=FuncMixNorm1Shatt(H,p)
    % where H is a LinOp object and p the order of the shatten norm (>=1)
    %
    % -- Properties
    % * |p|     order of the Shatten-norm (default 1)
    %
    % -- References
    % [1] Lefkimmiatis, S., Ward, J. P., & Unser, M. (2013). Hessian Schatten-norm regularization 
    %     for linear inverse problems. IEEE transactions on image processing, 22(5), 1873-1888.
    %
    % Please refer to the FUNC superclass for general documentation about
    % functional class
    % See also Func
	%
	%     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		p;       % order of the Shatten norm (>=1)
    end
    
    methods 
    	%% Constructor
        function this = FuncMixNorm1Shatt(H,p)
            this.name='Func MixNorm1-Shatten';
			this.isconvex= true; 
			if nargin==0 || isempty(H), 
				H=LinOpIdentity(); 
			else
				assert(length(H.sizeout)==3 && (H.sizeout(3)==3 || H.sizeout(3)==6),'sizeout of H should be [?,?,3 or 6]');
			end;
			if nargin<=1 || isempty(p), p=1; end;
			assert(p>=1,'p should be >=1');
			this.p=p
			this.set_H(H);
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
        	dim=size(x);
        	if dim(3)==3     % 2D
				[E,V]=svd2D_decomp(x);
			elseif dim(3)==6 % 3D
				[E,V]=svd3D_decomp(x);
			else
				error('third dimension of x should be 3 or 6');
			end
			if isinf(this.p)
				tmp=max(E,[],3);
			else
				tmp=sum(abs(E).^this.p,3).^(1/this.p);
			end
			y=sum(tmp(:));
        end
        %% Proximity operator of the functional
        function y=prox(this,x,alpha)
        	assert(isscalar(alpha),'alpha must be a scalar');
			% TODO IMPLEMENTS THE PROX (IF APPLICABLE)
        end
    end
end
