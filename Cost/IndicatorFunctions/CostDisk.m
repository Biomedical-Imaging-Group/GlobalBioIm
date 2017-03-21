classdef CostDisk < CostRing
    %% Circle indicator function
    %  Matlab Inverse Problems Library
    %
    % Obj = CostDisk(c):
    % Implement the proximity operator on the complex disk of radius c
    % $$ \phi(x) = 0 \textrm{ if } ||x||_2<=c  textrm{ and }  +\inf \textrm{ otherwise } $$
    % Following Matlab angle function if x=0 then CostDisk(x) = 1
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
        function this = CostDisk(radius)
            this.name='Cost Disk';
            
            if nargin >0
                assert(isnumeric(radius),'C should be numeric');
            else
                radius = 1;
            end
            
            if nargin==1
                sz = [];
            end 
            
            if isempty(sz)
                sz = size(radius);
            end
            set_H(this,sz,[]);
          
            
             this.outer = radius;
             this.inner = 0.;
        end
    end
end
