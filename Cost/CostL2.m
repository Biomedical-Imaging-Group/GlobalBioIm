classdef CostL2 < Cost
    %% CostL2 : Least Squares functional
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implement the cost function for the weighted L2 norm
    % $$ 1/2||Hx - y||^2_W $$
    % where H is a LinOp object (default LinOpIdentity), y are the data and W
    % is a weight matrix (LinOp, default LinOpIdentity).
    %
    % -- Example
    % F = CostL2(H,y,W);
    %
    % See also Cost, LinOp
    %
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch & Ferreol
    %     Soulez ferreol.soulez@epfl.ch
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
        W=1;          % weight matrix
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        fftHstardata=[]; % if LinOp is convolution, store the product conj(fftn(psf)).*fftn(data)
        wr;              % weighted residual
        isW=false;       % boolean true if a LinOp wght is given
    end
    
    methods
        %% Constructor
        function this = CostL2(H,y,wght)
            this.isconvex=true;
            % -- Set entries
            if nargin<2
                y=0;
            end
            if nargin<1
                H=[];
            end
            set_H(this,H,y);
            if nargin==3
                assert(isscalar(wght)||isa(wght,'LinOp'),'weight WGHT must be scalar or linop');
                this.W=wght;
                this.isW=true;
            end
            
            this.y=y;
            this.name='Cost L2';
            % -- Compute Lipschitz constant of the gradient (if the norm of H is known)
            if this.H.norm>=0;
                if isscalar(this.W)
                    this.lip=this.W.^2.*this.H.norm^2;
                else
                    if this.W.norm>=0
                        this.lip=this.H.norm^2*this.W.norm^2;
                    end
                end
            end
        end
        %% Evaluation of the Functional
        function f=eval(this,x)
            r=this.H.apply(x)-this.y;            
            this.wr=this.W*r;
            f=0.5*dot(r(:),this.wr(:));
        end
        %% Gradient of the Functional
        function g=grad(this,x)
            if nargin ==2
                this.wr = this.W*(this.H.apply(x)-this.y);
            end
            g = this.H.adjoint(this.wr) ;
        end
        %% Proximity operator of the functional
        function y=prox(this,x,alpha)
            assert(isscalar(alpha),'alpha must be a scalar');
            y=[];
            if isa(this.H,'LinOpIdentity')
                if this.isW && isa(this.W,'LinOpDiag')  % if weight is diagonal linop
                    y=(x+alpha*this.W.diag.*this.y)./(1+alpha.*this.W.diag);
                elseif ~this.isW % if no weight
                    y=(x+alpha*this.y)/(alpha+1);
                end
            elseif isa(this.H,'LinOpConv')  % if linop is convolution
                if isempty(this.fftHstardata)
                    this.fftHstardata=conj(this.H.mtf).*Sfft(this.y,this.H.Notindex);
                end
                if ~this.isW     % if no weight
                    y=iSfft((Sfft(x,this.H.Notindex) + alpha*this.fftHstardata)./(1+alpha*(abs(this.H.mtf).^2)), this.H.Notindex);
                    if ~this.H.iscomplex, y=real(y);end
                end
            end
            if isempty(y),error('Prox not implemented');end
        end
    end
end
