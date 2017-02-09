classdef LinOpNAME <  LinOp
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
        function this = LinOpNAME() %change the name of the constructor
            this.name ='';             
            this.issquare = false;     % is your operator square?
            this.iscomplex= true;      % is your operator complex?
            this.isinvertible = false; % true if the operator is invertible
            this.sizein = [];          % what is the size of the right hand side
            this.sizeout = [];         % what is the size of the left hand side
        end
        % MANDATORY METHODS
        function y = Apply(this,x)   
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
        end
%         FACULTATIVE METHODS
%         function y = Inverse(this,x)
%         end
%         function y = AdjointInverse(this,x)
%         end
%         function y = HtH(this,x) %  Apply the HtH matrix
%         end
%         function y = HHt(this,x) %  Apply the HHt matrix
%         end
%         function y = HtWH(this,x,W) %  Apply the HtWH matrix
%         end
    end
end
