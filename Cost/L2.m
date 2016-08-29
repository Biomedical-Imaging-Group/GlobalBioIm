classdef L2 < Cost
    %% L2 : L2 norm Cost Function 
    %  Matlab Inverse Problems Library
    % 
    % Obj = L2():
    % Implement the cost function for the weighted L2 norm
    % $$ \phi(x) = 1/2||x - a||^2_wght $$
    % 
   
    %  
    %     Copyright (C) 2016 F. Soulez ferreol.soulez@epfl.ch
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
        linOp   % 
        wght    % weight
        data    %
        x       %
        residual%
        cost    %
    end
    
    methods 
        function this = L2(linOp,data,wght)
            this.name='L2';
    end
    end
end
