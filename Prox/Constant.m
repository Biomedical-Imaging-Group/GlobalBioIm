classdef Constant < handle
    %% Constant : Constant Proximity Operator 
    %  Matlab Inverse Problems Library
    % 
    % Obj = Constant(cst):
    % Implement the proximity operator for the 
    % $$ \phi(x) = 0 \textrm{ if } x = cst textrm{ and }  +\inf \textrm{ otherwise } $$
    %
    %
    %% Properties
    % * |cst|         - constant (default 0)
    %
    %%
    
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
    
    properties (SetAccess = protected,GetAccess = public)
        cst
    end
    
    methods
        function this = Constant(cst)
            if nargin==0
                cst =0;
            end
             assert(isscalar(cst),'cst must be a scalar');
             this.cst = cst;
        end
        function y = Apply(this,x,~) % Apply the operator
            y = this.cst .* ones(size(x));
        end
        function y = Cost(this,x,~) % 
             y = 0.5 * norm(x - this.cst);
        end
        function y = FCost(~,~) 
            y = 0;
        end
        function y = Residual(this,x,rho) 
             y = (x - this.cst);
        end
    end
end
