classdef LinOpMatrix <  LinOp
    %% LinOpMatrix : Matrix operator
    %  Matlab Linear Operator Library 
    % Build the matrix operator such that Obj.apply(x) = M *x
    %
    % Example
    % Obj = LinOpMatrix(M)
    %
    % Be carefull it was not really tested
    %
    % Please refer to the LinOp superclass for documentation
    % See also LinOp
     
%     Copyright (C) 2015 F. Soulez  ferreol.soulez@epfl.ch
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
properties
    M
    index
    ndms
    sz
    rightsz
    leftsz
end
    methods
        function this = LinOpMatrix(M,index)
            this.name='Linop Matrix';
            this.M = M;
            this.sz =size(M);
            this.ndms = ndims(M);
            
                if isreal(M)
                    this.iscomplex= false;
                else
                    this.iscomplex= true;
                end
            
            if nargin == 1
                index = [];
            end
            if (~isempty(index))
                assert(isvector(index) && length(index)<= this.ndms && max(index)<= this.ndms,'The index should be a conformable  to sz');
                this.index = index;
            else
                 assert(this.ndms<3,'Without index parameter ndims(M) should be <3 ' );
                 this.index = 2;
            end
            
            T = true(this.ndms,1);
            T(this.index)=false;
            
            %size of the output = size of the input x length of the index
            % Special case for scalar vectors as matlab thought it is 2D matrix ;-(
            switch(length(this.index))
                case(this.ndms)
                this.sizeout= [1 1];
                this.sizein = this.sz;
                case(this.ndms-1)
                this.sizeout= [this.sz(T) 1];
                 this.sizein= [this.sz(~T) 1];
%                 this.leftsz= [ this.sz(T) 1];
%                 this.rightsz= [this.sz(~T) 1];
                otherwise
                this.sizeout= this.sz(T);
                this.sizein=  this.sz(~T);
            end
            this.rightsz = prod(this.sizein);
            this.leftsz = prod(this.sizeout);
            this.M= reshape(M,this.leftsz,this.rightsz);
            
        end
        function y = apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = reshape(this.M * reshape(x,[this.rightsz,1]),this.sizeout);
        end
        function y = adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
            y = reshape(this.M' * reshape(x,[this.leftsz,1]),this.sizein);
        end
    end
end

