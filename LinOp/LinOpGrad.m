classdef LinOpGrad <  LinOp
    %% LinOpGrad :  Finite difference operator
    %  Matlab Linear Operator Library
    %
    % -- Example
    % G = LinOpGrad(sz,index,bc,res)
    % Build the gradient operator (finite differences) to apply on a variable
    % of size SZ along the dimension indexed in INDEX (all by default)
    % The output is of size SZ x lenght(index)
    % bc corresponds to the boundary condition:
    %     - 'circular' (default)
    %     - 'zeros'
    %     - 'mirror'
    % res is a vector containing the resolution in each dimension (default: all 1)
    %      
    % NOTE: when circular boundary conditions are selected, the filter
    % corresponding to HtH (i.e. Laplacian) is available within the
    % attribute fHtH
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch, E. Soubies
    %     emmanuel.soubies@epfl.ch
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
        index;     % index along wich dimension are computed the finite differences
        lgthidx;   % length of INDEX
        ndms;      % number of dimension of the input
        bc;        % boundary condition (default mirror);
        res;       % resolution, vector of lenght ndms
        fHtH;      % Filter corresponding to HtH when using circular bc
    end
    methods
        function this = LinOpGrad(sz,index,bc,res)
            if nargin == 1
                index = [];
            end
            if nargin<=2 || isempty(bc)
                bc='circular';
            end
            if nargin<=3 || isempty(res)
            	res=ones(size(sz));
            end
            this.name ='LinOp Gradient';
            this.iscomplex= false;
            this.isinvertible=false;
            this.res=res;
            this.bc=bc;
            
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
            % size of the output = size of the input x length of the index
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizein(2) ==1
                this.sizeout= [this.sizein(1),this.lgthidx];
            else
                this.sizeout= this.sizein;
			end
			
			if this.lgthidx > 1
				this.sizeout(end+1) = this.lgthidx;
			end
            
            switch(this.ndms)
                case(1), this.norm=2/res(1);
                case(2), this.norm=2*sqrt(1/res(1)^2+1/res(2)^2);
                case(3), this.norm=2*sqrt(1/res(1)^2+1/res(2)^2+1/res(3)^2);
                case(4), this.norm=2*sqrt(1/res(1)^2+1/res(2)^2+1/res(3)^2+1/res(4)^2);
            end
            
            if strcmp(this.bc,'circular')
                % Set the filter for HtH
                this.fHtH=zeros(this.sizein);
                switch(this.ndms)
                    case(1), this.fHtH(1)=2;this.fHtH(2)=-1;this.fHtH(end)=-1;this.fHtH=this.fHtH/res(1)^2;
                    case(2), this.fHtH(1,1)=2/res(1)^2+2/res(2)^2;this.fHtH(1,2)=-1/res(2)^2;this.fHtH(2,1)=-1/res(1)^2;this.fHtH(1,end)=-1/res(2)^2;this.fHtH(end,1)=-1/res(1)^2;
                    case(3), this.fHtH(1,1,1)=2/res(1)^2+2/res(2)^2+2/res(3)^2;this.fHtH(1,2,1)=-1/res(2)^2;this.fHtH(2,1,1)=-1/res(1)^2;this.fHtH(1,end,1)=-1/res(2)^2;this.fHtH(end,1,1)=-1/res(1)^2;
                             this.fHtH(1,1,2)=-1/res(3)^2;this.fHtH(1,1,end)=-1/res(3)^2;
                    case(4), this.fHtH(1,1,1,1)=2/res(1)^2+2/res(2)^2+2/res(3)^2+2/res(4)^2;this.fHtH(1,2,1,1)=-1/res(2)^2;this.fHtH(2,1,1,1)=-1/res(1)^2;this.fHtH(1,end,1,1)=-1/res(2)^2;this.fHtH(end,1,1,1)=-1/res(1)^2;
                             this.fHtH(1,1,2,1)=-1/res(3)^2;this.fHtH(1,1,end,1)=-1/res(3)^2;this.fHtH(1,1,1,2)=-1/res(4)^2;this.fHtH(1,1,1,end)=-1/res(4)^2;
                end
            end
        end

        function y = apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = zeros(this.sizeout);
            nidx = 0;
            % switch according to the boundary condition
            switch(this.bc)
                case('mirror')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y(:,1) = (x([2:end,end])-x)/this.res(1);
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(:,:,nidx) = (x([2:end,end],:)-x)/this.res(1);
                                    case(2)
                                        y(:,:,nidx) = (x(:,[2:end,end])-x)/this.res(2);
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(:,:,:,nidx) = (x([2:end,end],:,:)-x)/this.res(1);
                                    case(2)
                                        y(:,:,:,nidx) = (x(:,[2:end,end],:)-x)/this.res(2);
                                    case(3)
                                        y(:,:,:,nidx) = (x(:,:,[2:end,end])-x)/this.res(3);
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(:,:,:,:,nidx) = (x([2:end,end],:,:,:)-x)/this.res(1);
                                    case(2)
                                        y(:,:,:,:,nidx) = (x(:,[2:end,end],:,:)-x)/this.res(2);
                                    case(3)
                                        y(:,:,:,:,nidx) = (x(:,:,[2:end,end],:)-x)/this.res(3);
                                    case(4)
                                        y(:,:,:,:,nidx) = (x(:,:,:,[2:end,end])-x)/this.res(4);
                                end
                            end
                    end
                case('circular')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y(:,1) = (x([2:end,1])-x)/this.res(1);
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(:,:,nidx) = (x([2:end,1],:)-x)/this.res(1);
                                    case(2)
                                        y(:,:,nidx) = (x(:,[2:end,1])-x)/this.res(2);
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(:,:,:,nidx) = (x([2:end,1],:,:)-x)/this.res(1);
                                    case(2)
                                        y(:,:,:,nidx) = (x(:,[2:end,1],:)-x)/this.res(2);
                                    case(3)
                                        y(:,:,:,nidx) = (x(:,:,[2:end,1])-x)/this.res(3);
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(:,:,:,:,nidx) = (x([2:end,1],:,:,:)-x)/this.res(1);
                                    case(2)
                                        y(:,:,:,:,nidx) = (x(:,[2:end,1],:,:)-x)/this.res(2);
                                    case(3)
                                        y(:,:,:,:,nidx) = (x(:,:,[2:end,1],:)-x)/this.res(3);
                                    case(4)
                                        y(:,:,:,:,nidx) = (x(:,:,:,[2:end,1])-x)/this.res(4);
                                end
                            end
                    end
                case('zeros')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y(1:end-1,1) = (x(2:end)-x(1:end-1))/this.res(1);
                            y(end,1) = -x(end)/this.res(1);
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1:end-1,:,nidx) = (x(2:end,:)-x(1:end-1,:))/this.res(1);
                                        y(end,:,nidx) = -x(end,:)/this.res(1);
                                    case(2)
                                        y(:,1:end-1,nidx) = (x(:,2:end)-x(:,1:end-1))/this.res(2);
                                        y(:,end,nidx) = -x(:,end)/this.res(2);
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1:end-1,:,:,nidx) = (x(2:end,:,:)-x(1:end-1,:,:))/this.res(1);
                                        y(end,:,:,nidx) = -x(end,:,:)/this.res(1);
                                    case(2)
                                        y(:,1:end-1,:,nidx) = (x(:,2:end,:)-x(:,1:end-1,:))/this.res(2);
                                        y(:,end,:,nidx) = -x(:,end,:)/this.res(2);
                                    case(3)
                                        y(:,:,1:end-1,nidx) = (x(:,:,2:end)-x(:,:,1:end-1))/this.res(3);
                                        y(:,:,end,nidx) = -x(:,:,end)/this.res(3);
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1:end-1,:,:,:,nidx) = (x(2:end,:,:,:)-x(1:end-1,:,:,:))/this.res(1);
                                        y(end,:,:,:,nidx) = -x(end,:,:,:)/this.res(1);
                                    case(2)
                                        y(:,1:end-1,:,:,nidx) = (x(:,2:end,:,:)-x(:,1:end-1,:,:))/this.res(2);
                                        y(:,end,:,:,nidx) = -x(:,end,:,:)/this.res(2);
                                    case(3)
                                        y(:,:,1:end-1,:,nidx) = (x(:,:,2:end,:)-x(:,:,1:end-1,:))/this.res(3);
                                        y(:,:,end,:,nidx) = -x(:,:,end,:)/this.res(3);
                                    case(4)
                                        y(:,:,:,1:end-1,nidx) = (x(:,:,:,2:end)-x(:,:,:,1:end-1))/this.res(4);
                                        y(:,:,:,end,nidx) = -x(:,:,:,end)/this.res(4);
                                end
                            end
                    end
            end
        end
        
        
        function y = adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            nidx = 0;
            y = zeros(this.sizein);
            switch(this.bc)
                case('mirror')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y= [[-x(1) ; (-x(2:end-1)+x(1:end-2))] ; x(end-1)]/this.res(1);
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:)       = y(1,:)  - x(1,:,nidx)/this.res(1);
                                        y(2:end-1,:) = y(2:end-1,:) + (-x(2:end-1,:,nidx)+x(1:end-2,:,nidx))/this.res(1);
                                        y(end,:)     = y(end,:) + x(end-1,:,nidx)/this.res(1);
                                    case(2)
                                        y(:,1)       = y(:,1)  - x(:,1,nidx)/this.res(2);
                                        y(:,2:end-1) = y(:,2:end-1) + (-x(:,2:end-1,nidx)+x(:,1:end-2,nidx))/this.res(2);
                                        y(:,end)     = y(:,end) + x(:,end-1,nidx)/this.res(2);
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:,:)       = y(1,:,:) - x(1,:,:,nidx)/this.res(1);
                                        y(2:end-1,:,:) = y(2:end-1,:,:) + (-x(2:end-1,:,:,nidx)+x(1:end-2,:,:,nidx))/this.res(1);
                                        y(end,:,:)     = y(end,:,:) + x(end-1,:,:,nidx)/this.res(1);
                                    case(2)
                                        y(:,1,:)       = y(:,1,:)  - x(:,1,:,nidx)/this.res(2);
                                        y(:,2:end-1,:) = y(:,2:end-1,:) + (-x(:,2:end-1,:,nidx)+x(:,1:end-2,:,nidx))/this.res(2);
                                        y(:,end,:)     = y(:,end,:) + x(:,end-1,:,nidx)/this.res(2);
                                    case(3)
                                        y(:,:,1)       = y(:,:,1)  - x(:,:,1,nidx)/this.res(3);
                                        y(:,:,2:end-1) = y(:,:,2:end-1) + (-x(:,:,2:end-1,nidx)+x(:,:,1:end-2,nidx))/this.res(3);
                                        y(:,:,end)     = y(:,:,end) + x(:,:,end-1,nidx)/this.res(3);
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:,:,:)       = y(1,:,:,:) - x(1,:,:,:,nidx)/this.res(1);
                                        y(2:end-1,:,:,:) = y(2:end-1,:,:,:) + (-x(2:end-1,:,:,:,nidx)+x(1:end-2,:,:,:,nidx))/this.res(1);
                                        y(end,:,:,:)     = y(end,:,:,:) + x(end-1,:,:,:,nidx)/this.res(1);
                                    case(2)
                                        y(:,1,:,:)       = y(:,1,:,:)  - x(:,1,:,:,nidx)/this.res(2);
                                        y(:,2:end-1,:,:) = y(:,2:end-1,:,:) + (-x(:,2:end-1,:,:,nidx)+x(:,1:end-2,:,:,nidx))/this.res(2);
                                        y(:,end,:,:)     = y(:,end,:,:) + x(:,end-1,:,:,nidx)/this.res(2);
                                    case(3)
                                        y(:,:,1,:)       = y(:,:,1,:)  - x(:,:,1,:,nidx)/this.res(3);
                                        y(:,:,2:end-1,:) = y(:,:,2:end-1,:) + (-x(:,:,2:end-1,:,nidx)+x(:,:,1:end-2,:,nidx))/this.res(3);
                                        y(:,:,end,:)     = y(:,:,end,:) + x(:,:,end-1,:,nidx)/this.res(3);
                                    case(4)
                                        y(:,:,:,1)       = y(:,:,:,1)  - x(:,:,:,1,nidx)/this.res(4);
                                        y(:,:,:,2:end-1) = y(:,:,:,2:end-1) + (-x(:,:,:,2:end-1,nidx)+x(:,:,:,1:end-2,nidx))/this.res(4);
                                        y(:,:,:,end)     = y(:,:,:,end) + x(:,:,:,end-1,nidx)/this.res(4);
                                end
                            end
                    end
                case('circular')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y= [[x(end)-x(1) ; (-x(2:end-1)+x(1:end-2))] ; (x(end-1)-x(end))]/this.res(1);
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:)       = y(1,:) + (x(end,:,nidx) - x(1,:,nidx))/this.res(1);
                                        y(2:end-1,:) = y(2:end-1,:) + (-x(2:end-1,:,nidx)+x(1:end-2,:,nidx))/this.res(1);
                                        y(end,:)     = y(end,:) + (x(end-1,:,nidx)-x(end,:,nidx))/this.res(1);
                                    case(2)
                                        y(:,1)       = y(:,1) + (x(:,end,nidx) - x(:,1,nidx))/this.res(2);
                                        y(:,2:end-1) = y(:,2:end-1) + (-x(:,2:end-1,nidx)+x(:,1:end-2,nidx))/this.res(2);
                                        y(:,end)     = y(:,end) + (x(:,end-1,nidx)-x(:,end,nidx))/this.res(2);
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:,:)       = y(1,:,:) + (x(end,:,:,nidx) - x(1,:,:,nidx))/this.res(1);
                                        y(2:end-1,:,:) = y(2:end-1,:,:) + (-x(2:end-1,:,:,nidx)+x(1:end-2,:,:,nidx))/this.res(1);
                                        y(end,:,:)     = y(end,:,:) + (x(end-1,:,:,nidx)-x(end,:,:,nidx))/this.res(1);
                                    case(2)
                                        y(:,1,:)       = y(:,1,:) + (x(:,end,:,nidx) - x(:,1,:,nidx))/this.res(2);
                                        y(:,2:end-1,:) = y(:,2:end-1,:) + (-x(:,2:end-1,:,nidx)+x(:,1:end-2,:,nidx))/this.res(2);
                                        y(:,end,:)     = y(:,end,:) + (x(:,end-1,:,nidx)-x(:,end,:,nidx))/this.res(2);
                                    case(3)
                                        y(:,:,1)       = y(:,:,1) + (x(:,:,end,nidx) - x(:,:,1,nidx))/this.res(3);
                                        y(:,:,2:end-1) = y(:,:,2:end-1) + (-x(:,:,2:end-1,nidx)+x(:,:,1:end-2,nidx))/this.res(3);
                                        y(:,:,end)     = y(:,:,end) + (x(:,:,end-1,nidx)-x(:,:,end,nidx))/this.res(3);
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:,:,:)       = y(1,:,:,:) +(x(end,:,:,:,nidx) - x(1,:,:,:,nidx))/this.res(1);
                                        y(2:end-1,:,:,:) = y(2:end-1,:,:,:) + (-x(2:end-1,:,:,:,nidx)+x(1:end-2,:,:,:,nidx))/this.res(1);
                                        y(end,:,:,:)     = y(end,:,:,:) + (x(end-1,:,:,:,nidx)-x(end,:,:,:,nidx))/this.res(1);
                                    case(2)
                                        y(:,1,:,:)       = y(:,1,:,:) +(x(:,end,:,:,nidx) - x(:,1,:,:,nidx))/this.res(2);
                                        y(:,2:end-1,:,:) = y(:,2:end-1,:,:) + (-x(:,2:end-1,:,:,nidx)+x(:,1:end-2,:,:,nidx))/this.res(2);
                                        y(:,end,:,:)     = y(:,end,:,:) + (x(:,end-1,:,:,nidx)-x(:,end,:,:,nidx))/this.res(2);
                                    case(3)
                                        y(:,:,1,:)       = y(:,:,1,:) +(x(:,:,end,:,nidx) - x(:,:,1,:,nidx))/this.res(3);
                                        y(:,:,2:end-1,:) = y(:,:,2:end-1,:) + (-x(:,:,2:end-1,:,nidx)+x(:,:,1:end-2,:,nidx))/this.res(3);
                                        y(:,:,end,:)     = y(:,:,end,:) + (x(:,:,end-1,:,nidx)-x(:,:,end,:,nidx))/this.res(3);
                                    case(4)
                                        y(:,:,:,1)       = y(:,:,:,1) +(x(:,:,:,end,nidx) - x(:,:,:,1,nidx))/this.res(4);
                                        y(:,:,:,2:end-1) = y(:,:,:,2:end-1) + (-x(:,:,:,2:end-1,nidx)+x(:,:,:,1:end-2,nidx))/this.res(4);
                                        y(:,:,:,end)     = y(:,:,:,end) + (x(:,:,:,end-1,nidx)-x(:,:,:,end,nidx))/this.res(4);
                                end
                            end
                    end
                case('zeros')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y= [[-x(1) ; (-x(2:end-1)+x(1:end-2))] ; x(end-1)-x(end)]/this.res(1);
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:)       = y(1,:)  - x(1,:,nidx)/this.res(1);
                                        y(2:end-1,:) = y(2:end-1,:) + (-x(2:end-1,:,nidx)+x(1:end-2,:,nidx))/this.res(1);
                                        y(end,:)     = y(end,:) + (x(end-1,:,nidx)-x(end,:,nidx))/this.res(1);
                                    case(2)
                                        y(:,1)       = y(:,1)  - x(:,1,nidx)/this.res(2);
                                        y(:,2:end-1) = y(:,2:end-1) + (-x(:,2:end-1,nidx)+x(:,1:end-2,nidx))/this.res(2);
                                        y(:,end)     = y(:,end) + (x(:,end-1,nidx)-x(:,end,nidx))/this.res(2);
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:,:)       = y(1,:,:)  - x(1,:,:,nidx)/this.res(1);
                                        y(2:end-1,:,:) = y(2:end-1,:,:) + (-x(2:end-1,:,:,nidx)+x(1:end-2,:,:,nidx))/this.res(1);
                                        y(end,:,:)     = y(end,:,:) + (x(end-1,:,:,nidx)-x(end,:,:,nidx))/this.res(1);
                                    case(2)
                                        y(:,1,:)       = y(:,1,:)  - x(:,1,:,nidx)/this.res(2);
                                        y(:,2:end-1,:) = y(:,2:end-1,:) + (-x(:,2:end-1,:,nidx)+x(:,1:end-2,:,nidx))/this.res(2);
                                        y(:,end,:)     = y(:,end,:) + (x(:,end-1,:,nidx)-x(:,end,:,nidx))/this.res(2);
                                    case(3)
                                        y(:,:,1)       = y(:,:,1)  - x(:,:,1,nidx)/this.res(3);
                                        y(:,:,2:end-1) = y(:,:,2:end-1) + (-x(:,:,2:end-1,nidx)+x(:,:,1:end-2,nidx))/this.res(3);
                                        y(:,:,end)     = y(:,:,end) + (x(:,:,end-1,nidx)-x(:,:,end,nidx))/this.res(3);
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                nidx = nidx +1;
                                switch(n)
                                    case(1)
                                        y(1,:,:,:)       = y(1,:,:,:)  - x(1,:,:,:,nidx)/this.res(1);
                                        y(2:end-1,:,:,:) = y(2:end-1,:,:,:) + (-x(2:end-1,:,:,:,nidx)+x(1:end-2,:,:,:,nidx))/this.res(1);
                                        y(end,:,:,:)     = y(end,:,:,:) + (x(end-1,:,:,:,nidx)-x(end,:,:,:,nidx))/this.res(1);
                                    case(2)
                                        y(:,1,:,:)       = y(:,1,:,:)  - x(:,1,:,:,nidx)/this.res(2);
                                        y(:,2:end-1,:,:) = y(:,2:end-1,:,:) + (-x(:,2:end-1,:,:,nidx)+x(:,1:end-2,:,:,nidx))/this.res(2);
                                        y(:,end,:,:)     = y(:,end,:,:) + (x(:,end-1,:,:,nidx)-x(:,end,:,:,nidx))/this.res(2);
                                    case(3)
                                        y(:,:,1,:)       = y(:,:,1,:)  - x(:,:,1,:,nidx)/this.res(3);
                                        y(:,:,2:end-1,:) = y(:,:,2:end-1,:) + (-x(:,:,2:end-1,:,nidx)+x(:,:,1:end-2,:,nidx))/this.res(3);
                                        y(:,:,end,:)     = y(:,:,end,:) + (x(:,:,end-1,:,nidx)-x(:,:,end,:,nidx))/this.res(3);
                                    case(4)
                                        y(:,:,:,1)       = y(:,:,:,1)  - x(:,:,:,1,nidx)/this.res(4);
                                        y(:,:,:,2:end-1) = y(:,:,:,2:end-1) + (-x(:,:,:,2:end-1,nidx)+x(:,:,:,1:end-2,nidx))/this.res(4);
                                        y(:,:,:,end)     = y(:,:,:,end) + (x(:,:,:,end-1,nidx)-x(:,:,:,end,nidx))/this.res(4);
                                end
                            end
                    end
            end
        end
        
        function y = HtH(this,x) %  apply the HtH matrix
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizein);
            y = zeros(this.sizein);
            switch(this.bc)
                case('mirror')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y = (2*x - x([1,1:end-1]) - x([2:end,end]))/this.res(1)^2;
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y = y + (2*x - x([1,1:end-1],:) - x([2:end,end],:))/this.res(1)^2;
                                    case(2)
                                        y = y + (2*x - x(:,[1,1:end-1]) - x(:,[2:end,end]))/this.res(2)^2;
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y = y + (2*x - x([1,1:end-1],:,:) - x([2:end,end],:,:))/this.res(1)^2;
                                    case(2)
                                        y = y + (2*x - x(:,[1,1:end-1],:) - x(:,[2:end,end],:))/this.res(2)^2;
                                    case(3)
                                        y = y + (2*x - x(:,:,[1,1:end-1]) - x(:,:,[2:end,end]))/this.res(3)^2;
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y = y + (2*x - x([1,1:end-1],:,:,:) - x([2:end,end],:,:,:))/this.res(1)^2;
                                    case(2)
                                        y = y + (2*x - x(:,[1,1:end-1],:,:) - x(:,[2:end,end],:,:))/this.res(2)^2;
                                    case(3)
                                        y = y + (2*x - x(:,:,[1,1:end-1],:) - x(:,:,[2:end,end],:))/this.res(3)^2;
                                    case(4)
                                        y = y + (2*x - x(:,:,:,[1,1:end-1]) - x(:,:,:,[2:end,end]))/this.res(4)^2;
                                end
                            end
                    end
                case('circular')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y = (2*x - x([end,1:end-1]) - x([2:end,1]))/this.res(1)^2;
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y = y + (2*x - x([end,1:end-1],:) - x([2:end,1],:))/this.res(1)^2;
                                    case(2)
                                        y = y + (2*x - x(:,[end,1:end-1]) - x(:,[2:end,1]))/this.res(2)^2;
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y = y + (2*x - x([end,1:end-1],:,:) - x([2:end,1],:,:))/this.res(1)^2;
                                    case(2)
                                        y = y + (2*x - x(:,[end,1:end-1],:) - x(:,[2:end,1],:))/this.res(2)^2;
                                    case(3)
                                        y = y + (2*x - x(:,:,[end,1:end-1]) - x(:,:,[2:end,1]))/this.res(3)^2;
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y = y + (2*x - x([end,1:end-1],:,:,:) - x([2:end,1],:,:,:))/this.res(1)^2;
                                    case(2)
                                        y = y + (2*x - x(:,[end,1:end-1],:,:) - x(:,[2:end,1],:,:))/this.res(2)^2;
                                    case(3)
                                        y = y + (2*x - x(:,:,[end,1:end-1],:) - x(:,:,[2:end,1],:))/this.res(3)^2;
                                    case(4)
                                        y = y + (2*x - x(:,:,:,[end,1:end-1]) - x(:,:,:,[2:end,1]))/this.res(4)^2;
                                end
                            end
                    end
                case('zeros')
                    % switch according to the number of dimension of the input
                    switch(this.ndms)
                        % 1 dimension
                        case(1)
                            y(1:end-1) = (2*x(1:end-1) - x([1,1:end-2]) - x(2:end))/this.res(1)^2;
                            y(end)=(2*x(end)-x(end-1))/this.res(1)^2;
                            % 2 dimensions
                        case(2)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y(1:end-1,:) = y(1:end-1,:) + (2*x(1:end-1,:) - x([1,1:end-2],:) - x(2:end,:))/this.res(1)^2;
                                        y(end,:)=y(end,:)+(2*x(end,:)-x(end-1,:))/this.res(1)^2;
                                    case(2)
                                        y(:,1:end-1) = y(:,1:end-1) + (2*x(:,1:end-1) - x(:,[1,1:end-2]) - x(:,2:end))/this.res(2)^2;
                                        y(:,end)=y(:,end)+(2*x(:,end)-x(:,end-1))/this.res(2)^2;
                                end
                            end
                            % 3 dimensions
                        case(3)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y(1:end-1,:,:) = y(1:end-1,:,:) + (2*x(1:end-1,:,:) - x([1,1:end-2],:,:) - x(2:end,:,:))/this.res(1)^2;
                                        y(end,:,:)=y(end,:,:)+(2*x(end,:,:)-x(end-1,:,:))/this.res(1)^2;
                                    case(2)
                                        y(:,1:end-1,:) = y(:,1:end-1,:) + (2*x(:,1:end-1,:) - x(:,[1,1:end-2],:) - x(:,2:end,:))/this.res(2)^2;
                                        y(:,end,:)=y(:,end,:)+(2*x(:,end,:)-x(:,end-1,:))/this.res(2)^2;
                                    case(3)
                                        y(:,:,1:end-1) = y(:,:,1:end-1) + (2*x(:,:,1:end-1) - x(:,:,[1,1:end-2]) - x(:,:,2:end))/this.res(3)^2;
                                        y(:,:,end)=y(:,:,end)+(2*x(:,:,end)-x(:,:,end-1))/this.res(3)^2;
                                end
                            end
                            % 4 dimensions
                        case(4)
                            for n=this.index
                                switch(n)
                                    case(1)
                                        y(1:end-1,:,:,:) = y(1:end-1,:,:,:) + (2*x(1:end-1,:,:,:) - x([1,1:end-2],:,:,:) - x(2:end,:,:,:))/this.res(1)^2;
                                        y(end,:,:,:)=y(end,:,:,:)+(2*x(end,:,:,:)-x(end-1,:,:,:))/this.res(1)^2;
                                    case(2)
                                        y(:,1:end-1,:,:) = y(:,1:end-1,:,:) + (2*x(:,1:end-1,:,:) - x(:,[1,1:end-2],:,:) - x(:,2:end,:,:))/this.res(2)^2;
                                        y(:,end,:,:)=y(:,end,:,:)+(2*x(:,end,:,:)-x(:,end-1,:,:))/this.res(2)^2;
                                    case(3)
                                        y(:,:,1:end-1,:) = y(:,:,1:end-1,:) + (2*x(:,:,1:end-1,:) - x(:,:,[1,1:end-2],:) - x(:,:,2:end,:))/this.res(3)^2;
                                        y(:,:,end,:)=y(:,:,end,:)+(2*x(:,:,end,:)-x(:,:,end-1,:))/this.res(3)^2;
                                    case(4)
                                        y(:,:,:,1:end-1) = y(:,:,:,1:end-1) + (2*x(:,:,:,1:end-1) - x(:,:,:,[1,1:end-2]) - x(:,:,:,2:end))/this.res(4)^2;
                                        y(:,:,:,end)=y(:,:,:,end)+(2*x(:,:,:,end)-x(:,:,:,end-1))/this.res(4)^2;
                                end
                            end
                    end
            end
            

        end
    end
    
end
