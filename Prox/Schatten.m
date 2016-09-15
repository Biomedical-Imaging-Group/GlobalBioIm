classdef Schatten < Prox
    %% Schatten Prox
    %  Matlab Inverse Problems Library
    % 
	% Shatten norm prox from Stamatis Lefkimmiatis (TIP, 2012)
    %
    %
    %% Properties
    %
    %%
    %     Copyright (C) 2015 michael.mccann@epfl.ch laurene.donati@epfl.ch
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
	  order; 
    end
    
    methods
        function this= Schatten(order)
            this.name='Schatten';
			this.order= order; 
        end
        function y = Apply(this,x,lambda) % Apply the operator
		  if ndims(x) == 3 % 2D case
			y = projectSpMat2x2(x,this.order,lambda);
		  elseif ndims(x) == 4 % 3D case
			y = projectSpMat3x3(x,this.order,lambda);
		  end 
		end
    end
end
