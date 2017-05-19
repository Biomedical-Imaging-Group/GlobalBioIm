classdef CostL1 < Cost
    %% CostL1 : CostL1 norm Proximity Operator
    %  Matlab Inverse Problems Library
    %
    % Obj = CostL1(H,y):
    % Implement the proximity operator for the L1 norm
    % $$ \phi(x) = |H.x-y| $$
    %
    % The option 'NonNegativity' add the non negativity constraint (default
    % false)
    
    %
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
    %%
    
    properties (SetAccess = protected,GetAccess = public)
        nonneg = false
    end
    
    methods
        function this = CostL1(H,y, varargin)
            this.name='Cost L1';
            this.isconvex=true;
            % -- Set entries
            if nargin<2
                y=[];
            end
            if nargin<1
                H=[];
            end
            set_H(this,H,y);
            
            for c=1:length(varargin)
                switch varargin{c}
                    case('NonNegativity')
                        this.nonneg = true;
                end
            end
            
        end
        function z = prox(this,x,alpha) % apply the operator
            assert(isscalar(alpha),'alpha must be a scalar');
            
            if this.H.isinvertible
                 tmp = this.H.apply(x)-this.y ;
                if this.nonneg
                    z =   this.H.inverse(max(tmp- alpha,0)+this.y);
                else
                    z =  this.H.inverse(sign(tmp) .* max( abs( tmp) - alpha,0)+this.y);
                end
            else
                error('Prox not implemented');
            end
        end
        function f=eval(this,x)
            if ~this.nonneg
                 tmp = this.H.apply(x)-this.y;
                 f=sum(abs(tmp(:)));
            else
                error('eval L1 not implemented');
            end
                 
        end       
    end
end
