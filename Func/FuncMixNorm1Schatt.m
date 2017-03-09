classdef FuncMixNorm1Schatt < Func
    %% FuncMixNorm1Schatt : Mixed l1-Schatten Norm
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
    % Note: The actual implementation works for Hx having one of the two following forms:
    %          - Hx (NxMx3) such that the Sp norm will be applied on each symetric 2x2
    %            matrix [Hx(n,m,1) Hx(n,m,2)
    %                    Hx(n,m,2), Hx(n,m,3)]
    %            and then the l1 norm on the two other dimensions
    %          - Hx (NxMxKx6) such that the Sp norm will be applied on each symetric 3x3
    %            matrix [Hx(n,m,k,1) Hx(n,m,k,2) Hx(n,m,k,3)
    %                    Hx(n,m,k,2) Hx(n,m,k,4)  Hx(n,m,k,5)
    %                    Hx(n,m,k,3) Hx(n,m,k,5) Hx(n,m,k,6)] 
    %            and then the l1 norm on the two other dimensions.
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
        function this = FuncMixNorm1Schatt(H,p)
            this.name='Func MixNorm1-Shatten';
			this.isconvex= true; 
			if nargin==0 || isempty(H), 
				H=LinOpIdentity(); 
			else
				assert((length(H.sizeout)==3 || length(H.sizeout)==4) && (H.sizeout(3)==3 || H.sizeout(4)==6),'sizeout of H should be [?,?,(?),3 or 6]');
			end;
			if nargin<=1 || isempty(p), p=1; end;
			assert(p>=1,'p should be >=1');
			this.p=p;
			this.set_H(H);
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
            dim=size(x);
            if dim(3)==3     % 2D
                [E,V]=svd2D_decomp(x);
            elseif dim(4)==6 % 3D
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
            dim=size(x);
            if this.p==1
                if dim(3)==3     % 2D
                    [E,V]=svd2D_decomp(x);
                    E=max(abs(E)-alpha,0).*sign(E);
                    y=svd2D_recomp(E,V);
                elseif dim(4)==6 % 3D
                    [E,V]=svd3D_decomp(x);
                    E=max(abs(E)-alpha,0).*sign(E);
                    y=svd3D_recomp(E,V);
                else
                    error('third dimension of x should be 3 or 6');
                end
            elseif this.p==2
                if dim(3)==3     % 2D
                    N=sqrt(x(:,:,1).^2+2*x(:,:,2).^2+x(:,:,3).^2);
                    y=repmat((N-1)./N,[1,1,3]).*x;
                elseif dim(4)==6 % 3D
                    N=sqrt(x(:,:,:,1).^2+2*x(:,:,:,2).^2+2*x(:,:,:,3).^2+x(:,:,:,4)+2*x(:,:,:,5)+x(:,:,:,6));
                    y=repmat((N-1)./N,[1,1,6]).*x;
                else
                    error('third dimension of x should be 3 or 6');
                end
            else
                error('prox not implemented');
            end
        end
    end
end
