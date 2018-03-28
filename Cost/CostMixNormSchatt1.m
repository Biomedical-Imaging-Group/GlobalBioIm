classdef CostMixNormSchatt1 < Cost
    % CostMixNormSchatt1 Mixed Schatten-l1 Norm [1]
    % $$C(\\mathrm{x}) :=   \\sum_n  \\| \\mathrm{x}_{n\\cdot} \\|_{Sp}, $$
    % for \\(p \\in [1,+\\infty]\\). Here, \\(\\|\\cdot\\|_{Sp}\\)  denotes the p-order Shatten norm 
    % defined by
    % $$ \\|\\mathrm{X}\\|_{Sp} = \\left[\\sum_k (\\sigma_k(\\mathrm{X}))^p\\right]^{1/p},$$
    % where \\(\\sigma_k(\\mathrm{X})\\) is the k-th singular value of \\(\\mathrm{X}\\). In other words it is the lp-norm
    % of the signular values of  \\(\\mathrm{X}\\).
    %
    % :param p: order of the Shatten norm (default 1)
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **Note** The actual implementation works for size (sz) having one of the two following forms: 
    %
    % - (NxMx3) such that the Sp norm will be applied on each symetric 2x2
    % $$ \\begin{bmatrix} \\mathrm{x}_{n m 1} & \\mathrm{x}_{n m 2} \\newline 
    % \\mathrm{x}_{n m 2} & \\mathrm{x}_{n m 3} \\end{bmatrix}$$
    % and then the \\(\\ell_1\\) norm on the two other dimensions.
    %
    % - (NxMxKx6) such that the Sp norm will be applied on each symetric 3x3
    % $$ \\begin{bmatrix} \\mathrm{x}_{n m k 1} & \\mathrm{x}_{n m k 2} & \\mathrm{x}_{n m k 3} \\newline 
    % \\mathrm{x}_{n m k 2} & \\mathrm{x}_{n m k 4} & \\mathrm{x}_{n m k 5} \\newline 
    % \\mathrm{x}_{n m k 3} & \\mathrm{x}_{n m k 5} & \\mathrm{x}_{n m k 6} \\newline  \\end{bmatrix}$$
    % and then the \\(\\ell_1\\) norm on the three other dimensions.
    %
    % **References**
    % [1] Lefkimmiatis, S., Ward, J. P., & Unser, M. (2013). Hessian Schatten-norm regularization
    % for linear inverse problems. IEEE transactions on image processing, 22(5), 1873-1888.
    %
    % **Example** C=CostMixNormSchatt1(sz,p,y)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`LinOp`
    
    %%    Copyright (C) 2017 
    %     E. Soubies emmanuel.soubies@epfl.ch
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
        reshDim; % index for reshaping
    end
    
    %% Constructor
    methods        
        function this = CostMixNormSchatt1(sz,p,y)
            % Verify if the mexgl files exist
            if (exist('svd2D_decomp')~=3)||(exist('svd2D_recomp')~=3)||(exist('svd3D_decomp')~=3)||(exist('svd3D_recomp')~=3)
                buildHessianSchatten();
            end
            
            if nargin<3, y=0; end
            this@Cost(sz,y);
            this.name='CostMixNormSchatt1';
            this.isConvex=true;
            this.isDifferentiable=false;
            if nargin<2 || isempty(p), p=1; end;
            assert(p>=1,'p should be >=1');
            assert(this.sizein(end)==3 || this.sizein(end)==6,'last dimension should be 3 or 6');
            this.p=p;   
            if this.sizein(end)==3
                this.reshDim=[this.sizein(1),prod(this.sizein(2:end-1)),3];
            elseif this.sizein(end)==6
                this.reshDim=[this.sizein(1),prod(this.sizein(2:end-1)),6];
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            dim=size(x);           
            if dim(end)==3     % 2D
                [E,~]=svd2D_decomp(reshape(x,this.reshDim));
            elseif dim(end)==6 % 3D
                [E,~]=svd3D_decomp(reshape(x,this.reshDim));
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
        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            dim=size(x);
            if this.p==1
                if dim(end)==3     % 2D
                    [E,V]=svd2D_decomp(reshape(x,this.reshDim));
                    E=max(abs(E)-alpha,0).*sign(E);
                    y=reshape(svd2D_recomp(E,V),this.sizein);
                elseif dim(end)==6 % 3D
                    [E,V]=svd3D_decomp(reshape(x,this.reshDim));
                    E=max(abs(E)-alpha,0).*sign(E);
                    y=reshape(svd3D_recomp(E,V),this.sizein);
                else
                    error('third dimension of x should be 3 or 6');
                end
            elseif this.p==2
                if dim(end)==3     % 2D
                    N=sqrt(x(:,:,1).^2+2*x(:,:,2).^2+x(:,:,3).^2);
                    y=repmat((N-1)./N,[ones(1,length(this.sizein)-1),3]).*x;
                elseif dim(end)==6 % 3D
                    N=sqrt(x(:,:,:,1).^2+2*x(:,:,:,2).^2+2*x(:,:,:,3).^2+x(:,:,:,4)+2*x(:,:,:,5)+x(:,:,:,6));
                    y=repmat((N-1)./N,[ones(1,length(this.sizein)-1),6]).*x;
                else
                    error('third dimension of x should be 3 or 6');
                end
            else
                y=applyProx_@Cost(this,x,alpha);
            end
        end
    end
end
