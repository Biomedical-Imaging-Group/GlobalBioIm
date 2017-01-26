classdef LinOpDiag <  LinOp
    %% LinOpDiag : Diagonal operator
    %  Matlab Linear Operator Library 
    %
    % Example:
    % Obj = LinOpDiag(diag)
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
        usecomplex = true;
    end
    methods
        function this = LinOpDiag(diag,varargin)
            this.name ='LinOp Diagonal';
            this.issquare = true;
            
          
            
            if isnumeric(diag)
                this.sizeout=size(diag);
                this.sizein=size(diag);
                if isreal(diag)
                    this.iscomplex= false;
                else
                    this.iscomplex= true;
                end
                
                  for c=1:length(varargin)
                switch varargin{c}
                    case('DontUseComplex')
                        this.usecomplex = false;
                         diag = reshape(diag,[],2);
                    diag = complex(diag(:,1),diag(:,2));
                end
                  end
            
                if all(diag)
                    this.isinvertible=true;
                else
                    this.isinvertible=false;
                end
                
                this.diag = diag;
            else
                error('diag value must be numeric');
            end           
            % -- Norm of the operator 
            this.norm=max(diag(:));
        end
        function y = Apply(this,x)
            if isequal(size(x),this.sizein)
                if ~this.usecomplex
                    x = reshape(x,[],2);
                    x = complex(x(:,1),x(:,2));
                    y =this.diag .* x;
                      y = cat(3,real(y),imag(y));
                    y = reshape(y , this.sizeout);
                else
                y =this.diag .* x;
                end
            else
                error('x should be the same size as diag: [%d, %d, %d, %d]',this.sizein);
            end
        end
        function y = Adjoint(this,x)
            
            if isequal(size(x),this.sizeout)
                if this.iscomplex
                    y =conj(this.diag) .*x;
                else
                if ~this.usecomplex
                    x = reshape(x,[],2);
                    x = complex(x(:,1),x(:,2));
                    y =conj(this.diag) .* x;
                      y = cat(3,real(y),imag(y));
                    y = reshape(y , this.sizeout);
                else
                    y =this.diag .*x;
                end
                end
            else
                error('x should be the same size as diag: [%d, %d, %d, %d]',this.sizein);
            end
        end
        
        function y = HtH(this,x) %  Apply the HtH matrix
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
