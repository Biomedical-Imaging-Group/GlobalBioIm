classdef TemplateCost < Cost
	% TODO: Put here the description of your Cost
    %
    % All attributes of parent class :class:`Cost` are inherited. 
	%
	% :param MyNewParam: TODO: details new parameter specific to this class
    %
    % See also :class:`Cost` :class:`LinOp`
    
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
        function this = TemplateCost(H,y)
        	% TODO SET THE INHERITED PROPERTIES (IF APPLICABLE)
            this.name='Cost NAME';
            this.isconvex= ???;
            this.lip= ???;
            % -- Set entries
            if nargin<2
                y=0;
            end
            if nargin<1
                H=[];
            end
            set_y(this,y);
            set_H(this,H);
			% TODO SET NEW DEFINED PROPERTIES
    	end

        function y=eval(this,x)
			% Reimplemented from parent class :class:`Cost`.
			
        end

		% TODO: Reimplements any method of class Cost that is relevant (prescise in method docstring any new specificities compared to class Cost. Otherwise put simply: Reimplemented from parent class :class:`Cost`.
		
    end
end
