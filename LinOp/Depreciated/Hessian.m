
classdef Hessian <  LinOp
    %% Hessian : 
    %  Matlab Linear Operator Library 
    % Build the identity operator such that x = Id.apply(x)
    %
    % Example
    % Id = Identity()
    %
    % Please refer to the LinOp superclass for documentation
    % See also LinOp
     
%     Copyright (C) 2015 F. Soulez  ferreol.soulez@epfl.ch
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
	  bc; %boundary conditions
    end

    methods
        function this = Hessian(sz, bc)
            this.name='Hessian';
            this.iscomplex=false;
            this.issquare=false;
			this.sizein=sz;
			this.bc=bc;
			if length(sz)==2 %2D case 
			  this.sizeout=[sz 3];
			elseif length(sz)==3 %3D case
			  this.sizeout=[sz 6];
			end
        end
        function y = apply(this,x)
		  if length(this.sizein)==2
			y = Hessian2D(x,this.bc);
		  elseif length(this.sizein)==3
			y = Hessian3D(x,this.bc);
		  end
        end
        function y = adjoint(this,x)
		  if length(this.sizein)==2 %2D case
			y = AdjHessian2D(x,this.bc);
		  elseif length(this.sizein)==3 %3D case
			y = AdjHessian3D(x,this.bc);
		  end
		end
    end
end

