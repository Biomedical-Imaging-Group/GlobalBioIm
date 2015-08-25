classdef Diagonal <  LinOp
    %% Diagonal : Diagonal operator
    %  Matlab Linear Operator Library 
    %
    % Example:
    % Obj = Diagonal(diag)
    %
    % Build the diagonal operator that multiply element wise the input by
    % the vector diag
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
        diag % diagonal vector
    end
    methods
        function this = Diagonal(diag)
            this.name ='Diagonal';
            this.issquare = true;
            if isnumeric(diag)
                this.diag = diag;
                this.sizeout=size(diag);
                this.sizein=size(diag);
                if isreal(diag)
                    this.iscomplex= false;
                else
                    this.iscomplex= true;
                end
                if all(diag)
                    this.isinvertible=true;
                else
                    this.isinvertible=false;
                end
            else
                error('diag value must be numeric');
            end
        end
        function y = Apply(this,x)
            if isequal(size(x),this.sizein)
                y =this.diag .* x;
            else
                error('x should be the same size as diag: [%d, %d, %d, %d]',this.sizein);
            end
        end
        function y = Adjoint(this,x)
            
            if isequal(size(x),this.sizeout)
                if this.iscomplex
                    y =conj(this.diag) .*x;
                else
                    y =this.diag .*x;
                end
            else
                error('x should be the same size as diag: [%d, %d, %d, %d]',this.sizein);
            end
        end
        function y = Inverse(this,x)
            if ( ~this.isinvertible)
                error('Operator non invertible');
            end
            y =(1./this.diag) .*x;
        end
        function y = AdjointInverse(this,x)
            if ( ~this.isinvertible)
                error('Operator non invertible');
            end
            if this.iscomplex
                y =conj(1./this.diag) .*x;
            else
                y = (1./this.diag) .*x;
            end
        end
    end
end
