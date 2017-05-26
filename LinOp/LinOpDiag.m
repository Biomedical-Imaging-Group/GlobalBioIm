classdef LinOpDiag <  LinOp
    %% LinOpDiag : Diagonal operator
    %  Matlab Linear Operator Library
    %
    % Example:
    % Obj = LinOpDiag(diag)
    %
    % Build the diagonal operator that multiplies element wise the input by
    % the vector DIAG  or
    %
    % Obj = LinOpDiag(diag,sz)
    %
    % Build the diagonal operator that multiplies element wise an input of size SZ by
    % the scalar DIAG
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
        function this = LinOpDiag(diag,sz)
            this.name ='LinOp Diagonal';            
            if ~isnumeric(diag)
                error('diag must be numeric');
            end
            
            % collapse repeated diagonal element to a scalar
            % Very slow
            %if length(unique(diag)) == 1
            %    diag = unique(diag);
            %end
            
            if isscalar(diag)
                if isempty(sz) || ~issize(sz)
                    error('must provide sz argument when diag is a scalar');
                end
                this.sizein = sz;
                this.sizeout = sz;
            else
                this.sizeout=size(diag);
                this.sizein=size(diag);
            end
            
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
            
            this.diag = diag;
            
            % -- Norm of the operator
            this.norm=max(abs(diag(:)));
        end
        function y = apply(this,x)
            if isequal(size(x),this.sizein)
                y =this.diag .* x;
            else
                error('x should be the same size as diag: [%d, %d, %d, %d]',this.sizein);
            end
        end
        function y = adjoint(this,x)
            
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
        
        function y = HtH(this,x) %  apply the HtH matrix
            if isequal(size(x),this.sizeout)
                if this.iscomplex
                    y =abs(this.diag).^2 .*x;
                else
                    y =this.diag.^2 .*x;
                end
            else
                error('x should be the same size as diag: [%d, %d, %d, %d]',this.sizein);
            end
        end
        function y=HHt(this,x)
            y=this.HtH(x);
        end
        
        function y = inverse(this,x)
            if ( ~this.isinvertible)
                error('Operator non invertible');
            end
            y =(1./this.diag) .*x;
        end
        function y = adjointInverse(this,x)
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
