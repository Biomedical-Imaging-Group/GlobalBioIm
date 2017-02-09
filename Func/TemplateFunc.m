classdef FuncNAME < Func
    %% FuncNAME : TODO ADD DESCRIPTION
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % TODO ADD A DESCRIPTION
    %
    % -- Example
    % TODO ADD INSTANTIATION EXAMPLE
    %
    % -- Properties
    % TODO ADD NEW PROPERTIES
    %
    % -- Methods
    % TODO ADD NEW METHODS
    %
    % -- References
    % TODO ADD REFERENCES IF NEEDED
    %
    % Please refer to the FUNC superclass for general documentation about
    % functional class
    % See also Func
	%
    %     Copyright (C) TODO YEAR NAME EMAIL
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
		%% TODO SET HERE NEW PROTECTED SET AND PUBLIC READ PROPERTIES
		%% IF NEEDED.
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
		%% TODO SET HERE NEW FULLY PROTECTED PROPERTIES 
		%% (E.G. INTERNAL VARIABLE USED TO AVOID MULTIPLE COMPUTATION)
    end
    
    methods 
    	%% Constructor
        function this = FuncNAME(~)
        	% TODO SET THE INHERITED PROPERTIES (IF APPLICABLE)
            this.name='Func NAME';
			this.isconvex= ???; 
			this.lip= ???;
			% call this.set_H (will also set the sizein parameter)
			% TODO SET NEW DEFINED PROPERTIES
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			% TODO IMPLEMENT THE EVAL METHOD
        end
        %% Gradient of the Functional
        function g=grad(this,x)
			% TODO IMPLEMENTS THE GRADIENT (IF APPLICABLE)
        end
        %% Proximity operator of the functional
        function y=prox(this,x,alpha)
        	assert(isscalar(alpha),'alpha must be a scalar');
			% TODO IMPLEMENTS THE PROX (IF APPLICABLE)
        end
        %% Function that set properly the operator H (has to be modified if new properties is???H are added)
        function set_H(this,H)
        	% TODO REIMPLEMENTS THE METHOD SET_H (IF NEEDED)
        end
    end
end
