classdef L2 < handle
    %% L2 : L2 norm Proximity Operator 
    %  Matlab Inverse Problems Library
    % 
    % Obj = L2(a):
    % Implement the proximity operator for the L2 norm
    % $$ \phi(x) = ||x - a||^2 $$
    % 
   
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
        a
    end
    
    methods 
        function this = L2(a)
            if isempty(a)
                a =0;
            end
            assert(isnumeric(a),'alpha must be a numeric');
            if isscalar(a)
                this.a = a;
            else
                this.a = a;
                this.sz = size(a);
            end
        end
        function y = Apply(this,x,alpha) % Apply the operator
            assert(isscalar(alpha),'alpha must be a scalar');
            this.alpha = alpha;
            if ~isscalar(this.a)
               assert( isequal(size(x),this.size),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein); 
            end
            y = 1./(2.*this.alpha + 1) * (2.*this.alpha.* x - this.a);
        end
    end
end
