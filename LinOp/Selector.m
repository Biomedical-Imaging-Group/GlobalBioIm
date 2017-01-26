classdef Selector <  LinOp
    %% Selector : Selector operator
    %  Matlab Linear Operator Library
    %
    % Example:
    % Obj = Selector(sel)
    %
    % Build the operator selecting values of the input vector according to
    % the boolean value of sel
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
        sel % diagonal vector
    end
    methods
        function this = Selector(sel,varargin)
            p = inputParser;
            addOptional(p,'KeepDimensions',false); %true or false.
            parse(p,varargin{:});
            
            this.name ='Selector';
            this.issquare = false;
            this.iscomplex= true;
            assert(islogical(sel),'The input selector should be boolean');
            
            if true(p.Results.KeepDimensions)
                   y = max(sum(sel,1));         %
                   x =max(sum(sel,2));
                   this.sizeout=[sum(y), sum(x)];
            else
                    this.sizeout=[sum(sel(:)), 1];
                    this.sel = sel;
            end
            this.sizein=size(sel);
            this.isinvertible=false;
        end
        
        function y = Apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y =x(this.sel);
            y = reshape(y,this.sizeout);
        end
        
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            y = zeros(this.sizein);
            y(this.sel) = x;
        end
        function y = HtH(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = x;
            y(~this.sel) = 0;
        end
        function y = HHt(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            y = x;
        end
    end
end
