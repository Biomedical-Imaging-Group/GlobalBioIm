classdef CostKullLeib < Cost
    %% CostKullLeib : KullbackLeibler divergence
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implements the Kullback-Leibler divergence between H.x and y .
    % $$ \sum_n -y_n log((H*x)_n + bet) + (H*x)_n,
    % where H is a LinOp object (default LinOpIdentity), y are the data
    % and bet is a scalar (default 0) to smooth the function at zero.
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
        Hx;     % H.apply(x)
    end
    
    methods
        %% Constructor
        function this = CostKullLeib(H,y,bet)
            this.name='Cost Kullback-Leibler';
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
                this.bet=bet;
            end
            tmp=ones(this.H.sizeout);
            this.He1=this.H.adjoint(tmp);
            
            % -- Compute Lipschitz constant of the gradient (if the norm of H is known)
            if ((this.bet>0)&&(this.H.norm>=0));
                this.lip=max(this.y(:))./this.bet^2*this.H.norm^2;
            end
        end
        %% Evaluation of the Functional
        function f=eval(this,x)
            this.Hx=this.H.apply(x);
            if ~any(this.Hx(:)<0)
                f=Inf;
            else
                if (this.bet~=0)
                    f=sum(-this.y(:).*log(this.Hx(:)+this.bet) + this.Hx(:));
                else
                    ft = zeros(size(this.Hx));
                    zidx = (this.Hx~=0);
                    ft(zidx)=-this.y(zidx).*log(this.Hx(zidx)) + this.Hx(zidx);
                    f=sum(ft(:));
                end
            end
        end
        %% Gradient of the Functional
        function g=grad(this,x)
            if nargin ==2
            this.Hx=this.H.apply(x);
            end
            g= this.He1 - this.H.adjoint(this.y./(this.Hx+this.bet));
        end
        %% Proximity operator of the functional
        function z=prox(this,x,alpha)
            z=[];
            if isa(this.H,'LinOpIdentity') % if operator identity
                if (this.bet~=0)
                    delta=(x-alpha-this.bet).^2+4*(x*this.bet + alpha*(this.y-this.bet));
                    z=zeros(size(x));
                    mask=delta>=0;
                    z(mask)=0.5*(x(mask)-alpha-this.bet + sqrt(delta(mask)));
                end
            end
            if isempty(z),error('Prox not implemented');end
        end
    end
end
