classdef LinOpSelector <  LinOp
   %% LinOpSelector : Selector operator
    %  Matlab Linear Operator Library
    %
    % Example:
    % Obj = LinOpSelector(sel) or Obj = Selector(sel,'KeepDimensions',true)
    %
    % Builds the operator by selecting the values of the input vector according to
    % the boolean value of sel. If the option KeepDimensions is activated,
    % sel has to select a compact rectangle. The compact rectangle is not stored for memory
    % reasons.
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
        sel         % diagonal vector, boolean
        XLIM
        YLIM
        KeepDimensions %true or false
    end
    methods
        function this = LinOpSelector(sel,varargin)
            p = inputParser;
            addOptional(p,'KeepDimensions',false); %true or false.
            parse(p,varargin{:});

             this.KeepDimensions = p.Results.KeepDimensions;
            this.name ='LinOp Selector';
            this.issquare = false;
            this.iscomplex= true;
            assert(islogical(sel),'The input selector should be boolean');
            
			this.norm = 1;
			
            if true(p.Results.KeepDimensions)
                   [row,col]= find(sel ~=0);
                   this.XLIM = [min(col) max(col)];          %x index, start and end positions of the selector
                   this.YLIM =  [min(row) max(row)];      %y index, start and end positions of the selector
                   this.sizeout=[this.XLIM(2)-this.XLIM(1)+1, this.YLIM(2)-this.YLIM(1)+1];
            else
                    this.sizeout=[sum(sel(:)), 1];
                    this.sel = sel;
            end
            this.sizein=size(sel);
            this.isinvertible=false;
        end
        
        function y = apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            if true(this.KeepDimensions)
                y = x(this.YLIM(1):this.YLIM(2),this.XLIM(1):this.XLIM(2));
            else
                y =x(this.sel);
            end
            y = reshape(y,this.sizeout);
        end
        
        function y = adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            y = zeros(this.sizein);
            if true(this.KeepDimensions)
                y(this.YLIM(1):this.YLIM(2),this.XLIM(1):this.XLIM(2)) = x;
            else
                y(this.sel) = x;
            end
        end
        function y = HtH(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
                if true(this.KeepDimensions)
                    y = zeros(this.sizein);
                    y(this.YLIM(1):this.YLIM(2),this.XLIM(1):this.XLIM(2)) = x(this.YLIM(1):this.YLIM(2),this.XLIM(1):this.XLIM(2));
                else
                    y = x;
                    y(~this.sel) = 0;                
                end
        end
        function y = HHt(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            y = x;
        end
    end
end
