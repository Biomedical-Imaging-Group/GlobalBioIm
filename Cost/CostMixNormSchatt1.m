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
    %   * (NxMx3) such that the Sp norm will be applied on each symetric 2x2
    %     $$ \\begin{bmatrix} \\mathrm{x}_{n m 1} & \\mathrm{x}_{n m 2} \\newline 
    %     \\mathrm{x}_{n m 2} & \\mathrm{x}_{n m 3} \\end{bmatrix}$$
    %     and then the \\(\\ell_1\\) norm on the two other dimensions.
    %   * (NxMxKx6) such that the Sp norm will be applied on each symetric 3x3
    %     $$ \\begin{bmatrix} \\mathrm{x}_{n m k 1} & \\mathrm{x}_{n m k 2} & \\mathrm{x}_{n m k 3} \\newline 
    %     \\mathrm{x}_{n m k 2} & \\mathrm{x}_{n m k 4} & \\mathrm{x}_{n m k 5} \\newline 
    %     \\mathrm{x}_{n m k 3} & \\mathrm{x}_{n m k 5} & \\mathrm{x}_{n m k 6} \\newline  \\end{bmatrix}$$
    %     and then the \\(\\ell_1\\) norm on the three other dimensions.
    %
    % **References**
    % [1] Lefkimmiatis, S., Ward, J. P., & Unser, M. (2013). Hessian Schatten-norm regularization
    % for linear inverse problems. IEEE transactions on image processing, 22(5), 1873-1888.
    %
    % **Example** C=CostMixNormSchatt1(sz,p,y)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`LinOp`
    
    %% GUI-Header
    % GUInotation-mnSch1-
    % GUIcall-CostMixNormSchatt1(InputSize,order,data)-
    % GUIparam-InputSize-vecInt-[]-Input size of the cost function (e.g. [512 512]). If empty, the size of the data will be used.
    % GUIparam-data-file-[]-data vector (default 0)
    % GUIparam-order-vecInt-[]-order of the Shatten norm (default 1) 
    
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
                this.reshDim=[this.sizein(1),prod(this.sizein(2:end-1)),1,6];
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            [E,~]=this.svdDecomp(x);
            if isinf(this.p)
                tmp=max(E,[],3);
            else
                tmp=sum(abs(E).^this.p,3).^(1/this.p);
            end
            y=sum(tmp(:));
        end
        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            
            if this.p==1
                [E,V]=this.svdDecomp(x);
                E=max(abs(E)-alpha,0).*sign(E);
                y=reshape(this.svdRecomp(E,V),this.sizein);
            elseif this.p==2
                szx=size(x);
                dim=floor(sqrt(2*szx(end)));
                diag_inds=cumsum(dim+1:-1:2)-dim;
                x2=x.^2;
                Frob=sqrt(2*sum(x2,3) - sum(x2(:,:,diag_inds),3));
                y=max(1-alpha./Frob,0).*x;
            else
                y=applyProx_@Cost(this,x,alpha);
            end
        end
    end
    
    
    %% Internal methods
    methods (Access = protected)
        function [E,V]=svdDecomp(this,x)   
            
            global isGPU
            dim=size(x);      % 2D
            if dim(end)==3
                if isGPU==1
                    tt=(abs(x(:,:,2))<eps);
                    E=zeros_([dim(1:2),2]);
                    trace=x(:,:,1)+x(:,:,3);
                    delta=(x(:,:,1)-x(:,:,3)).^2+4*x(:,:,2).^2;
                    E(:,:,1)=0.5*(trace+sqrt(delta));
                    E(:,:,2)=0.5*(trace-sqrt(delta));
                    V=zeros_([dim(1:2),2]);
                    n=sqrt((E(:,:,1)-x(:,:,1)).^2+x(:,:,2).^2)+tt;
                    V(:,:,1)=x(:,:,2)./n.*(1-tt) + tt;
                    V(:,:,2)=(E(:,:,1)-x(:,:,1))./n.*(1-tt);
                else
                    [E,V]=svd2D_decomp(reshape(x,this.reshDim));
                    %[E,V]=svd2D_decomp(x);
                end
            elseif dim(end)==6 % 3D
                if isGPU==1
                    error([this.name,' cannot be used with GpuArray (used mex files are not supported by GpuArray). Instead, you can use CudaMat']);
                end
                [E,V]=svd3D_decomp(reshape(x,this.reshDim));
              %  [E,V]=svd3D_decomp(x);
            else
                error('last dimension of x should be 3 or 6');
            end
        end
        
        function x=svdRecomp(this,E,V)  
            
            global isGPU
            dim=size(E);
            if dim(end)==2     % 2D
                if isGPU==1
                    x=zeros_([dim(1:2),3]);
                    x(:,:,1)=E(:,:,1).*V(:,:,1).^2 + E(:,:,2).*V(:,:,2).^2;
                    x(:,:,2)=V(:,:,1).*V(:,:,2).*(E(:,:,1)-E(:,:,2));
                    x(:,:,3)=E(:,:,1).*V(:,:,2).^2+E(:,:,2).*V(:,:,1).^2;
                else
                    x=svd2D_recomp(E,V);
                end
            elseif dim(end)==3 % 3D
                if isGPU==1
                    error([this.name,' cannot be used with GpuArray (used mex files are not supported by GpuArray). Instead, you can use CudaMat']);
                end
                x=svd3D_recomp(E,V);
            else
                error('third dimension of E and V should be 3 or 6');
            end
        end
    end
end
