classdef CostLeastSquares < Cost
    %% CostLeastSquares : Least Squares functional
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implement the cost function for the weighted L2 norm
    % $$ 1/2||Hx - d||^2_W $$
    % where H is a LinOp object (default LinOpIdentity), d are the data and W
    % is a weight matrix (LinOp, default LinOpIdentity).
    %
    % -- Example
    % F = CostLeastSquares(d,[],W);
    %
    % Please refer to the FUNC superclass for general documentation about
    % functional class
    % See also Cost, LinOp
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
        data           % data
        W=[];          % weight matrix
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        WplusWt=[];      % sum of W plus its transpose
        Hd               % application of the adjoint to data
        fftHstardata=[]; % if LinOp is convolution, store the product conj(fftn(psf)).*fftn(data)
        isW=false;       % boolean true if a LinOp wght is given
    end
    
    methods
        %% Constructor
        function this = CostLeastSquares(data,H,wght)
            this.isconvex=true;
            % -- Set entries
            if nargin==1 || isempty(H)
                H=LinOpIdentity(size(data));
            end
            if nargin==3
                this.W=wght;
                this.WplusWt=wght+wght';
                this.isW=true;
            end
            this.data=data;
            this.name='Cost Least Squares';
            this.set_H(H);
            assert( isequal(size(data),this.H.sizeout),'H sizeout and data size are not equal');
            % -- Compute Lipschitz constant of the gradient (if the norm of H is known)
            if this.H.norm>=0;
                if this.isW
                    if this.W.norm>=0
                        this.lip=this.H.norm^2*this.W.norm^2;
                    end
                else
                    this.lip=this.H.norm^2;
                end
            end
        end
        %% Evaluation of the Functional
        function y=eval(this,x)
            r=this.H.Apply(x)-this.data;
            if this.isW
                wr=this.W.Apply(r);
                y=0.5*dot(r(:),wr(:));
            else
                y=0.5*norm(r(:))^2;
            end
        end
        %% Gradient of the Functional
        function g=grad(this,x)
            if this.isW
                g = 0.5*(this.H.HtWH(x,this.WplusWt)-this.Hd);
            else
                g = this.H.HtH(x) - this.Hd;
            end
        end
        %% Proximity operator of the functional
        function y=prox(this,x,alpha)
            assert(isscalar(alpha),'alpha must be a scalar');
            y=[];
            if isa(this.H,'LinOpIdentity')
                if this.isW && isa(this.W,'LinOpDiag')  % if weight is diagonal linop
                    y=(x+alpha*this.W.diag.*this.data)./(1+alpha.*this.W.diag);
                elseif ~this.isW % if no weight
                    y=(x+alpha*this.data)/(alpha+1);
                end
            elseif isa(this.H,'LinOpConv')  % if linop is convolution
                if isempty(this.fftHstardata)
                    this.fftHstardata=conj(this.H.mtf).*Sfft(this.data,this.H.Notindex);
                    if ~this.H.iscomplex, this.fftHstardata=real(this.fftHstardata);end
                end
                if ~this.isW     % if no weight
                    y=iSfft((Sfft(x,this.H.Notindex) + alpha*this.fftHstardata)./(1+alpha*(abs(this.H.mtf).^2)), this.H.Notindex);
                    if ~this.H.iscomplex, y=real(y);end
                end
            end
            if isempty(y),error('Prox not implemented');end
        end
        %% Function that set properly the operator H (has to be modified if new properties is???H are added)
        function set_H(this,H)
            this.H=H;
            this.sizein=this.H.sizein;
            if this.isW
                this.Hd=this.H.Adjoint(this.WplusWt.Apply(this.data));
            else
                this.Hd=this.H.Adjoint(this.data);
            end
        end
    end
end
