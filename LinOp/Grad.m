classdef Grad <  LinOp
    %% Grad :  Finite difference operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = Grad(sz,index)
    % Finite operator operator
    % Compute finite differences of a vector of size SZ along the dimension
    % indexed in INDEX (all by default)
    % the output is zero padded to have size conformable with the input
    % The output is then of size SZ x lenght(index)
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
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
        index % index along wich dimension are computed the finite differences
        lgthidx % length of INDEX
        ndms  % number of dimension of the input
    end
    methods
        function this = Grad(sz,index)
            if nargin == 1
                index = [];
            end
            this.name ='Grad';
            this.iscomplex= false;
            this.isinvertible=false;
            
            
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizein = sz;
            this.ndms = length(this.sizein);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizein(2) ==1
                this.ndms = 1;
            end
            
            if (~isempty(index))
                assert(isvector(index) && length(index)<= this.ndms && max(index)<= this.ndms,'The index should be a conformable  to sz');
                this.index = index;
            else
                this.index = 1:this.ndms;
            end
            this.lgthidx = length(this.index);
            %size of the output = size of the input x length of the index
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizein(2) ==1
                this.sizeout= [this.sizein(1),this.lgthidx];
            else
                this.sizeout= [this.sizein,this.lgthidx];
            end
        end
        function y = Apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = zeros(this.sizeout);
            nidx = 0;
            % switch according to the number of dimension of the input
            switch(this.ndms)
                % 1 dimension
                case(1)
                    % y(:,1) = padarray(x(1:end-1,1) - x(2:end,1),1,'post');
                    y(:,1) = padarray(diff(x,1,1),1,'post');
                    % 2 dimensions
                case(2)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y(:,:,nidx) = padarray(diff(x,1,1),1,'post');
                            case(2)
                                y(:,:,nidx) = padarray(diff(x,1,2),[0 1],'post');
                        end
                    end
                    % 3 dimensions
                case(3)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y(:,:,:,nidx) = padarray(diff(x,1,1),1,'post');
                            case(2)
                                y(:,:,:,nidx) = padarray(diff(x,1,2),[0 1],'post');
                            case(3)
                                y(:,:,:,nidx) = padarray(diff(x,1,3),[0 0 1],'post');
                        end
                    end
                    % 4 dimensions
                case(4)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y(:,:,:,:,nidx) =  padarray(diff(x,1,1),1,'post');
                            case(2)
                                y(:,:,:,:,nidx) = padarray(diff(x,1,2),[0 1],'post');
                            case(3)
                                y(:,:,:,:,nidx) = padarray(diff(x,1,3),[0 0 1],'post');
                            case(4)
                                y(:,:,:,:,nidx) = padarray(diff(x,1,4),[0 0  0 1],'post');
                        end
                    end
            end
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            nidx = 0;
            y = zeros(this.sizein);
            % switch according to the number of dimension of the input
            switch(this.ndms)
                % 1 dimension
                case(1)
                    y =  padarray( x(1:end-1),1,'pre') -  padarray(  x(1:end-1),1,'post');
                    % 2 dimensions
                case(2)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y = y +  padarray( x(1:end-1,:,nidx),1,'pre') - padarray(  x(1:end-1,:,nidx),1,'post');
                            case(2)
                                y = y +  padarray( x(:,1:end-1,nidx),[0 1],'pre')- padarray(  x(:,1:end-1,nidx),[0 1],'post');
                        end
                    end
                    % 3 dimensions
                case(3)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y = y - padarray(  x(1:end-1,:,:,nidx),1,'post') + padarray( x(1:end-1,:,:,nidx),1,'pre');
                            case(2)
                                y = y- padarray(  x(:,1:end-1,:,nidx),[0 1],'post') + padarray( x(:,1:end-1,:,nidx),[0 1],'pre');
                            case(3)
                                y = y - padarray(  x(:,:,1:end-1,nidx),[0  0 1],'post') + padarray( x(:,:,1:end-1,nidx),[0 0 1],'pre');
                        end
                    end
                    % 4 dimensions
                case(4)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y = y - padarray(  x(1:end-1,:,:,:,nidx),1,'post') + padarray( x(1:end-1,:,:,:,nidx),1,'pre');
                            case(2)
                                y = y - padarray(  x(:,1:end-1,:,:,nidx),[0 1],'post') + padarray( x(:,1:end-1,:,:,nidx),[0 1],'pre');
                            case(3)
                                y = y - padarray(  x(:,:,1:end-1,:,nidx),[0  0 1],'post') + padarray( x(:,:,1:end-1,:,nidx),[0 0 1],'pre');
                            case(4)
                                y = y - padarray(  x(:,:,:,1:end-1,nidx),[0  0 0 1],'post') + padarray( x(:,:,:,1:end-1,nidx),[0 0 0 1],'pre');
                        end
                    end
            end
            
            
        end
        function y = Gram(this,x) %  Apply the Gram matrix
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizein);
            nidx = 0;
            y = zeros(this.sizein);
            % switch according to the number of dimension of the input
            switch(this.ndms)
                % 1 dimension
                case(1)
                    y = - padarray( x(1:end-1),1,-x(1),'pre') -  padarray(  x(2:end),1,-x(end),'post') + 2*padarray(  x(2:end-1),1,'both');
                    % 2 dimensions
                case(2)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y = y + padarray( x(1:end-1,:),1,'post') -  padarray(  x(2:end,:),1,'post')  - padarray( x(1:end-1,:),1,'pre') +  padarray(  x(2:end,:),1,'pre');
                            case(2)
                                y = y + padarray( x(:,1:end-1),[0 1],'post') -  padarray(  x(:,2:end),[0 1],'post')  - padarray( x(:,1:end-1),[0 1],'pre') +  padarray(  x(:,2:end),[0 1],'pre');
                        end
                    end
                    % 3 dimensions
                case(3)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y = y + padarray( x(1:end-1,:,:),1,'post') -  padarray(  x(2:end,:,:),1,'post')  - padarray( x(1:end-1,:,:),1,'pre') +  padarray(  x(2:end,:,:),1,'pre');
                            case(2)
                                y = y + padarray( x(:,1:end-1,:),[0 1],'post') -  padarray(  x(:,2:end,:),[0 1],'post')  - padarray( x(:,1:end-1,:),[0 1],'pre') +  padarray(  x(:,2:end,:),[0 1],'pre');
                            case(3)
                                y = y + padarray( x(:,:,1:end-1),[0 0 1],'post') -  padarray(  x(:,:,2:end),[0 0 1],'post')  - padarray( x(:,:,1:end-1),[0 0 1],'pre') +  padarray(  x(:,:,2:end),[0 0 1],'pre');
                        end
                    end
                    % 4 dimensions
                case(4)
                    for n=this.index
                        nidx = nidx +1;
                        switch(n)
                            case(1)
                                y = y + padarray( x(1:end-1,:,:,:),1,'post') -  padarray(  x(2:end,:,:,:),1,'post')  - padarray( x(1:end-1,:,:,:),1,'pre') +  padarray(  x(2:end,:,:,:),1,'pre');
                            case(2)
                                y = y + padarray( x(:,1:end-1,:,:),[0 1],'post') -  padarray(  x(:,2:end,:,:),[0 1],'post')  - padarray( x(:,1:end-1,:,:),[0 1],'pre') +  padarray(  x(:,2:end,:,:),[0 1],'pre');
                            case(3)
                                y = y + padarray( x(:,:,1:end-1,:),[0 0 1],'post') -  padarray(  x(:,:,2:end,:),[0 0 1],'post')  - padarray( x(:,:,1:end-1,:),[0 0 1],'pre') +  padarray(  x(:,:,2:end,:),[0 0 1],'pre');
                            case(4)
                                y = y + padarray( x(:,:,:,1:end-1),[0 0 0 1],'post') -  padarray(  x(:,:,:,2:end),[0 0 0 1],'post')  - padarray( x(:,:,:,1:end-1),[0 0 0 1],'pre') +  padarray(  x(:,:,:,2:end),[0 0 0 1],'pre');
                                
                        end
                    end
            end
            
            
            
            
            
        end
    end
    
end
