classdef CostCircle < CostRing
    %% Circle indicator function
    %  Matlab Inverse Problems Library
    %
    % Obj = CostCircle(c):
    % Implement the proximity operator on the complex circle of radius c
    % $$ \phi(x) = 0 \textrm{ if } ||x||_2=c  textrm{ and }  +\inf \textrm{ otherwise } $$
    % Following Matlab angle function if x=0 then Circle(x) = 1
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
        function this = CostCircle(radius)
            this.name='Cost Circle';
            
            set_H(this,[],[]);
          
            
            if nargin >0
                if isscalar(radius)
                    this.radius = radius;
                else
                    if isnumeric(radius)
                        this.radius = radius;
                        this.sizein = size(radius);
                    else
                        error('C should be numeric');
                    end
                end
            else
                radius = 1;
            end
             this.outer = radius;
             this.inner = radius;
        end
    end
end
