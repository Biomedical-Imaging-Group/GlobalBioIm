classdef FuncMultScalar < Func
    %% FuncMultScalar : Multiply a Func by a scalar
    %  Matlab Inverse Problems Library
    %
    % -- Example
    % F = FuncSum(Func,alpha)
    % Multiply the Func by the scalar alpha
    %
    % Please refer to the FUNC superclass for general documentation about
    % functional class
    % See also Func
	%
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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
 
    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
        func;      % Func
        alpha;     % scalar factor
    end
    
    methods 
    	%% Constructor
        function this = FuncMultScalar(func,alpha)
            this.name='Func Multiply Scalar';
			this.func = func;
			this.alpha=alpha;
			this.sizein =  this.func.sizein;
			this.isconvex=func.isconvex; 
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			y=this.alpha*this.func.eval(x);
        end
        %% Gradient of the Functional
        function g=grad(this,x)
			g=this.alpha*this.func.grad(x);
        end
    end
end
