classdef NonNegativity < handle
    %% NonNegativity : NonNegative Proximity Operator 
    %  Matlab Inverse Problems Library
    % 
    % Obj = NonNegativity():
    % Implement the proximity operator for the non negativity constraint
    % $$ \phi(x) = 0 \textrm{ if } x \ge 0 textrm{ and }  +\inf \textrm{ otherwise } $$
    %
    %
    %% Properties
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
    end
    
    methods
        function y = Apply(~,x,~) % Apply the operator
            y = max(x,0.);
        end
        function Cost(~,x,~) % 
            y = norm(max(x,0.));
        end
        function y = FCost(~,~,~) 
            y = 0;
        end
        function y = Residual(~,x,~) 
             y = min(x,0.);
        end
    end
end
