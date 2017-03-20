classdef CostKullLeib < Cost
    %% CostKullLeib : KullbackLeibler divergence
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implements the Kullback-Leibler divergence between H.x and y .
    % $$ \sum_n -y_n log((H*x)_n + bet) + (H*x)_n,
    % where H is a LinOp object (default LinOpIdentity), y are the data
    % and bet is a small scalar (default 1e-10) to smooth the function at zero.
    %
    % -- Example
    % F=CostKullLeib(H,y,bet)
    %
    % Please refer to the COST superclass for general documentation about
    % functional class
    % See also Cost
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
        bet= 0;     % smoothing parameter, if bet=0 then the unsmoothed version is used
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        He1;     % Application of the adjoint of H to a vector of ones
    end
    
    methods
        %% Constructor
        function this = CostKullLeib(H,y,bet)
            this.name='Cost Kullback-Leibler';
            this.isconvex=true;
            % -- Set entries
            if nargin>0
                if ~isempty(H)
                    this.set_H(H);
                end
            end
            if nargin>1
                if ~isempty(y)
                    this.y=y;
                end
            end
            
            if nargin==3
                this.bet=bet;
            end
            if ~isempty(this.H.sizeout)
                assert((~isscalar(this.y)) && isequal(size(this.y),this.H.sizeout),'H sizeout and y size are not equal');
                tmp=ones(this.sizein);
                this.He1=this.H.Adjoint(tmp);
            end
            % -- Compute Lipschitz constant of the gradient (if the norm of H is known)
            if ((this.bet>0)&&(this.H.norm>=0));
                this.lip=max(this.y(:))./this.bet^2*this.H.norm^2;
            end
        end
        %% Evaluation of the Functional
        function f=eval(this,x)
            if isempty(this.H.sizeout)
                tmp=ones(this.sizein);
                this.He1=this.H.Adjoint(tmp);
            end
            tmp=this.H.Apply(x);
            assert(any(tmp(:)<0) ,'H.x must be non-negative');
            if (all(this.bet(:)))
                f=sum(-this.y(:).*log(tmp(:)+this.bet) + tmp(:));
            else
                f = zeros(size(tmp));
                zidx = (tmp(:)~=0);
                f(zidx)=sum(-this.y(zidx).*log(tmp(zidx)+this.bet(zidx)) + tmp(zidx));
            end
            
        end
        %% Gradient of the Functional
        function g=grad(this,x)
            if isempty(this.H.sizeout)
                tmp=ones(this.sizein);
                this.He1=this.H.Adjoint(tmp);
            end
            g= this.He1 - this.H.Adjoint(this.y./(this.H.Apply(x)+this.bet));
        end
        %% Proximity operator of the functional
        function z=prox(this,x,alpha)
            z=[];
            if isa(this.H,'LinOpIdentity') % if operator identity
                delta=(x-alpha-this.bet).^2+4*(x*this.bet + alpha*(this.y-this.bet));
                z=zeros(size(x));
                mask=delta>=0;
                z(mask)=0.5*(x(mask)-alpha-this.bet + sqrt(delta(mask)));
            end
            if isempty(z),error('Prox not implemented');end
        end
    end
end
