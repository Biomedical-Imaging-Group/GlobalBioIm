classdef CostL2 < Cost
    % Weighted L2 norm cost function
    % $$C(\\mathrm{x}) := \\frac12\\|\\mathrm{Hx} - \\mathrm{y}\\|^2_W $$
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param W: weighting :class:`LinOpDiag` object or scalar (default :class:`LinOpIdentity`)
    %
    % See also :class:`Cost` :class:`LinOp`

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
        W=1;    % weight matrix      
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        fftHstardata=[]; % if LinOp is convolution, store the product conj(fftn(psf)).*fftn(data)
        isW=false;       % boolean true if a LinOp wght is given
    end
    
    methods
        function this = CostL2(H,y,wght)
            this.isconvex=true;
            % -- Set entries
            if nargin<2
                y=0;
            end
            if nargin<1
                H=[];
            end
            set_y(this,y);
            set_H(this,H);
            
            if nargin==3
                assert(isscalar(wght)||isa(wght,'LinOpDiag'),'weight WGHT must be scalar or Diagonal LinOp');
                this.W=wght;
                this.isW=true;
            end
            
            this.y=y;
            this.name='Cost L2';
            % -- Compute Lipschitz constant of the gradient (if the norm of H is known)
            if this.H.norm>=0;
                if isnumeric(this.W)
                    this.lip=this.W.*this.H.norm^2;
                else
                    if this.W.norm>=0
                        this.lip=this.H.norm^2*this.W.norm;
                    end
                end
            end
        end
        
        function f=eval(this,x)
        	% Reimplemented from parent class :class:`Cost`.
        	
            if(isscalar(this.y)&&(this.y==0))
                r=this.H.apply(x);
            else
                r=this.H.apply(x)-this.y;
            end
            wr=this.W*r;
            f=0.5*dot(r(:),wr(:));
        end

        function g=grad(this,x)
        	% Reimplemented from parent class :class:`Cost`.
        	% $$ \\nabla C(\\mathrm{x}) = \\mathrm{H^* W (Hx - y)} $$
        	% It is L-Lipschitz continuous with \\( L \\leq \\|\\mathrm{H}\\|^2 \\|\\mathrm{W}\\|\\).
        	
            if(isscalar(this.y)&&(this.y==0))
                r=this.H.apply(x);
            else
                r=this.H.apply(x)-this.y;
            end
            wr=this.W*r;
            g = this.H.adjoint(wr) ;
        end
        
        function [cost , gradient] = eval_grad(this,x)
        	% Reimplemented from parent class :class:`Cost`.
            
            if(isscalar(this.y)&&(this.y==0))
                r=this.H.apply(x);
            else
                r=this.H.apply(x)-this.y;
            end
            wr=this.W*r;
            cost=0.5*dot(r(:),wr(:));
            gradient = this.H.adjoint(wr);
        end

        function y=prox(this,x,alpha)
        	% Reimplemented from parent class :class:`Cost` if
        	%
        	% - the operator :attr:`H`  is a :class:`LinOpIdentity`,
        	% $$\\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\frac{\\mathrm{x}+\\alpha \\mathrm{W}\\mathrm{y}}{1+\\alpha \\mathrm{W}}$$
        	% where the division is component-wise.
        	%
        	% - the operator :attr:`H`  is a :class:`LinOpConv` and :attr:`W`  is a :class:`LinOpIdentity`;
        	% $$\\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathcal{F}^{-1}\\left(\\frac{\\mathcal{F}(\\mathrm{x}) + \\alpha  \\mathcal{F}(\\mathrm{H}^*)\\mathcal{F}(\\mathrm{y})  }{1+\\alpha \\vert\\mathcal{F}(\\mathrm{H})\\vert^2} \\right)$$
        	% where \\(\\mathcal{F} \\) stands for the Fourier transform.
        	
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
            elseif isa(this.H,'LinOpSelector')
                if ~this.isW
                    t=this.H.apply(x);
                    y=x + this.H.adjoint((t+alpha*this.y)/(alpha+1) - t);
                else
                end
            end
            if isempty(y),error('Prox not implemented');end
        end
    end
end
