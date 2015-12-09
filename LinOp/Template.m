classdef Template <  LinOp
    %% Template : Template for Linear  Operator
    %  Matlab Linear Operator Library 
    %
    % TODO Write your do here
    % Example:
    %
    % Please refer to the LinOp superclass for documentation
    % See also LinOp
    
    
    
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
       % Add your extra properties
    end
    methods
        function this = Template() %change the name of the constructor
            this.name ='';             %
            this.issquare = true;      % is your operator square?
            this.iscomplex= true;      % is your operator complex?
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
