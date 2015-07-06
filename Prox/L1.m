classdef L1 < Prox
    %% L1 : L1 norm Proximity Operator 
    %  Matlab Inverse Problems Library
    % 
    % Obj = L1():
    % Implement the proximity operator for the L1 norm
    % $$ \phi(x) = |x| $$
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
        function this = L1(varargin)
            this.name='L1';
            for c=1:length(varargin)
                switch varargin{c}
                    case('NonNegativity')
                        this.nonneg = true;
                end
            end
            
        end
        function y = Apply(this,x,alpha) % Apply the operator
            assert(isscalar(alpha),'alpha must be a scalar');
            this.alpha = alpha;
            if this.nonneg
                y =  max( x - alpha,0);
            else
                y = sign(x) .* max( abs(x) - alpha,0);
            end
        end
    end
end
