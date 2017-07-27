classdef TemplateLinOp <  LinOp
    %% LinOpNAME : Template for Linear  Operator
    %  Matlab Linear Operator Library 
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
    % Please refer to the LinOp superclass for documentation
    % See also LinOp   
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

    properties (SetAccess = protected,GetAccess = public)
       % Add your extra properties
    end
    methods
        function this = TemplateLinOp() %change the name of the constructor
            this.name ='';             
            this.issquare = false;     % is your operator square?
            this.isComplex= true;      % is your operator complex?
            this.isInvertible = false; % true if the operator is invertible
            this.sizein = [];          % what is the size of the right hand side
            this.sizeout = [];         % what is the size of the left hand side
		end
	end
	
	methods (Access = protected)
        % MANDATORY METHODS
        function y = apply_(this,x)   
			y = [];
		end
		
        function y = adjoint_(this,x)
           y = [];
		end
		
%         OPTIONAL METHODS
%         function y = inverse_(this,x)
%         end
%         function y = adjointInverse_(this,x)
%         end
%         function y = HtH_(this,x) %  apply the HtH matrix
%         end
%         function y = HHt_(this,x) %  apply the HHt matrix
%         end
%         function y = HtWH_(this,x,W) %  apply the HtWH matrix
%         end
    end
end
