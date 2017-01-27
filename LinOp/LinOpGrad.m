classdef LinOpGrad <  LinOp
    %% LinOpGrad :  Finite difference operator
    %  Matlab Linear Operator Library
    %
    % -- Example
    % G = LinOpGrad(sz,index,res)
    % Finite operator operator
    % Compute finite differences of a vector of size SZ along the dimension
    % indexed in INDEX (all by default)
    % the output is zero padded to have size conformable with the input
    % The output is then of size SZ x lenght(index)
    % res is a vector containing the resolution in each dimension (default: all 1)
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
        index;     % index along wich dimension are computed the finite differences
        lgthidx;   % length of INDEX
        ndms;      % number of dimension of the input
        res;       % resolution, vector of lenght ndms
    end
    methods
        function this = LinOpGrad(sz,index,res)
            if nargin == 1
                index = [];
            end
            if nargin==2 || isempty(res)
            	res=ones(size(sz));
            end
            this.name ='LinOp Gradient';
            this.iscomplex= false;
            this.isinvertible=false;
            this.issquare = false;
            this.res=res;
            
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
				case(4),	this.norm=2*sqrt(1/res(1)^2+1/res(2)^2+1/res(3)^2+1/res(4)^2);		
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
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            nidx = 0;
            y = zeros(this.sizein);
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
                            	y(1,:)       = y(1,:) - x(1,:,nidx)/this.res(1);
                            	y(2:end-1,:) = y(2:end-1,:) + (-x(2:end-1,:,nidx)+x(1:end-2,:,nidx))/this.res(1);
                            	y(end,:)     = y(end,:) + x(end-1,:,nidx)/this.res(1);
                            case(2)
                            	y(:,1)       = y(:,1) - x(:,1,nidx)/this.res(2);
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
                            	y(:,1,:)       = y(:,1,:) - x(:,1,:,nidx)/this.res(2);
                            	y(:,2:end-1,:) = y(:,2:end-1,:) + (-x(:,2:end-1,:,nidx)+x(:,1:end-2,:,nidx))/this.res(2);
                            	y(:,end,:)     = y(:,end,:) + x(:,end-1,:,nidx)/this.res(2);
                            case(3)
                            	y(:,:,1)       = y(:,:,1) - x(:,:,1,nidx)/this.res(3);
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
                            	y(:,1,:,:)       = y(:,1,:,:) - x(:,1,:,:,nidx)/this.res(2);
                            	y(:,2:end-1,:,:) = y(:,2:end-1,:,:) + (-x(:,2:end-1,:,:,nidx)+x(:,1:end-2,:,:,nidx))/this.res(2);
                            	y(:,end,:,:)     = y(:,end,:,:) + x(:,end-1,:,:,nidx)/this.res(2);
                            case(3)
                            	y(:,:,1,:)       = y(:,:,1,:) - x(:,:,1,:,nidx)/this.res(3);
                            	y(:,:,2:end-1,:) = y(:,:,2:end-1,:) + (-x(:,:,2:end-1,:,nidx)+x(:,:,1:end-2,:,nidx))/this.res(3);
                            	y(:,:,end,:)     = y(:,:,end,:) + x(:,:,end-1,:,nidx)/this.res(3); 
                            case(4)
                            	y(:,:,:,1)       = y(:,:,:,1) - x(:,:,:,1,nidx)/this.res(4);
                            	y(:,:,:,2:end-1) = y(:,:,:,2:end-1) + (-x(:,:,:,2:end-1,nidx)+x(:,:,:,1:end-2,nidx))/this.res(4);
                            	y(:,:,:,end)     = y(:,:,:,end) + x(:,:,:,end-1,nidx)/this.res(4);       
                        end
                    end
            end
            
            
        end
        function y = HtH(this,x) %  Apply the HtH matrix
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizein);
            nidx = 0;
            y = zeros(this.sizein);
            % switch according to the number of dimension of the input
            switch(this.ndms)
                % 1 dimension
                case(1)
                	y = (2*x - x([1,1:end-1]) - x([2:end,end]))/this.res(1)^2;
                % 2 dimensions
                case(2)
                    for n=this.index
                        nidx = nidx +1;
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
                        nidx = nidx +1;
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
                        nidx = nidx +1;
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

        end
    end
    
end
