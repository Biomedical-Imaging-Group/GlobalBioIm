classdef LinOpHess <  LinOp
    %% LinOpHess : Hessian operator 
    %  Matlab Linear Operator Library 
    %
    % -- Example
    % Hess = LinOpHess(sz)
    % Build the Hessian operator to apply on a variable os size SZ
    %
    % Note: Only 2D and 3D cases are implemented.
    %       The output of the Apply() method is:
    %         - for 2D case: a [sz,3] matrix containing [d^2F/dxx;d^2F/dxy;d^2F/dyy] 
    %		  - for 3D case: a [sz,6] matrix containing [d^2F/dxx;d^2F/dxy;d^2F/dxz;d^2F/dyy;d^2F/dyz;d^2F/dzz]
    %       These size are the input sizes for the Adjoint() method
    %
    % 
    % Please refer to the LinOp superclass for documentation
    % See also LinOp
    % 
	%     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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
	  ndms;      % number of dimension of the input
    end

    methods
        function this = LinOpHess(sz)
            this.name='LinOp Hessian';
            this.iscomplex=false;
            this.issquare=false;
			this.sizein=sz;
			this.ndms = length(this.sizein);
			switch(this.ndms)
				case(2), this.sizeout=[sz 3];
				case(3), this.sizeout=[sz 6];
			end
        end
        function y = Apply(this,x)
			assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = zeros(this.sizeout);
            nidx = 0;
            % switch according to the number of dimension of the input
            switch(this.ndms)
                % 2 dimension
                case(2)
                	% xx
					y(1:end-2,:,1)=x(3:end,:)-2*x(2:end-1,:)+x(1:end-2,:);
					y(end-1,:,1)=x(end-1,:)-2*x(end,:);
					y(end,:,1)=x(end,:);
					% xy
					y(1:end-1,1:end-1,2)=x(2:end,2:end) - x(2:end,1:end-1)-x(1:end-1,2:end) + x(1:end-1,1:end-1);
					y(end,1:end-1,2)=x(end,1:end-1)-x(end,2:end);
					y(1:end-1,end,2)=x(1:end-1,end)-x(2:end,end);
					y(end,end,2)=x(end,end);
					% yy
					y(:,1:end-2,3)=x(:,3:end)-2*x(:,2:end-1)+x(:,1:end-2);
					y(:,end-1,3)=x(:,end-1)-2*x(:,end);
					y(:,end,3)=x(:,end);
					%}
                % 3 dimensions
                case(3)
                    % xx
					y(1:end-2,:,:,1)=x(3:end,:,:)-2*x(2:end-1,:,:)+x(1:end-2,:,:);
					y(end-1,:,:,1)=x(end-1,:,:)-2*x(end,:,:);
					y(end,:,:,1)=x(end,:,:);
					% xy
					y(1:end-1,1:end-1,:,2)=x(2:end,2:end,:) - x(2:end,1:end-1,:)-x(1:end-1,2:end,:) + x(1:end-1,1:end-1,:);
					y(end,1:end-1,:,2)=x(end,1:end-1,:)-x(end,2:end,:);
					y(1:end-1,end,:,2)=x(1:end-1,end,:)-x(2:end,end,:);
					y(end,end,:,2)=x(end,end,:);
					% xz
					y(1:end-1,:,1:end-1,3)=x(2:end,:,2:end) - x(2:end,:,1:end-1)-x(1:end-1,:,2:end) + x(1:end-1,:,1:end-1);
					y(end,:,1:end-1,3)=x(end,:,1:end-1)-x(end,:,2:end);
					y(1:end-1,:,end,3)=x(1:end-1,:,end)-x(2:end,:,end);
					y(end,:,end,3)=x(end,:,end);
					% yy
					y(:,1:end-2,:,4)=x(:,3:end,:)-2*x(:,2:end-1,:)+x(:,1:end-2,:);
					y(:,end-1,:,4)=x(:,end-1,:)-2*x(:,end,:);
					y(:,end,:,4)=x(:,end,:);
					% yz
					y(:,1:end-1,1:end-1,5)=x(:,2:end,2:end) - x(:,2:end,1:end-1)-x(:,1:end-1,2:end) + x(:,1:end-1,1:end-1);
					y(:,end,1:end-1,5)=x(:,end,1:end-1)-x(:,end,2:end);
					y(:,1:end-1,end,5)=x(:,1:end-1,end)-x(:,2:end,end);
					y(:,end,end,5)=x(:,end,end);
					% zz
					y(:,:,1:end-2,6)=x(:,:,3:end)-2*x(:,:,2:end-1)+x(:,:,1:end-2);
					y(:,:,end-1,6)=x(:,:,end-1)-2*x(:,:,end);
					y(:,:,end,6)=x(:,:,end);
            end
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            nidx = 0;
            y = zeros(this.sizein);
            % switch according to the number of dimension of the input
            switch(this.ndms)
                % 2 dimension
                case(2)
                	% xx
					y(1,:)=y(1,:) + x(1,:,1);
					y(2,:)=y(2,:) + x(2,:,1) - 2*x(1,:,1);
					y(3:end,:)=y(3:end,:) + x(1:end-2,:,1) - 2*x(2:end-1,:,1) + x(3:end,:,1);
					% xy
					y(2:end,2:end)=y(2:end,2:end) + x(1:end-1,1:end-1,2) - x(1:end-1,2:end,2) - x(2:end,1:end-1,2) + x(2:end,2:end,2);
					y(2:end,1) = y(2:end,1) -x(1:end-1,1,2)+ x(2:end,1,2);
					y(1,2:end) = y(1,2:end) - x(1,1:end-1,2) + x(1,2:end,2);
					y(1,1)=y(1,1)+x(1,1,2);
					% yy
					y(:,1)=y(:,1) + x(:,1,3);
					y(:,2)=y(:,2) + x(:,2,3) - 2*x(:,1,3);
					y(:,3:end)=y(:,3:end) + x(:,1:end-2,3) - 2*x(:,2:end-1,3) + x(:,3:end,3);
					%}
                % 3 dimensions
                case(3)
                    % xx
					y(1,:,:)=y(1,:,:) + x(1,:,:,1);
					y(2,:,:)=y(2,:,:) + x(2,:,:,1) - 2*x(1,:,:,1);
					y(3:end,:,:)=y(3:end,:,:) + x(1:end-2,:,:,1) - 2*x(2:end-1,:,:,1) + x(3:end,:,:,1);
					% xy
					y(2:end,2:end,:)=y(2:end,2:end,:) + x(1:end-1,1:end-1,:,2) - x(1:end-1,2:end,:,2) - x(2:end,1:end-1,:,2) + x(2:end,2:end,:,2);
					y(2:end,1,:) = y(2:end,1,:) -x(1:end-1,1,:,2)+ x(2:end,1,:,2);
					y(1,2:end,:) = y(1,2:end,:) - x(1,1:end-1,:,2) + x(1,2:end,:,2);
					y(1,1,:)=y(1,1,:)+x(1,1,:,2);
					% xz
					y(2:end,:,2:end)=y(2:end,:,2:end) + x(1:end-1,:,1:end-1,3) - x(1:end-1,:,2:end,3) - x(2:end,:,1:end-1,3) + x(2:end,:,2:end,3);
					y(2:end,:,1) = y(2:end,:,1) -x(1:end-1,:,1,3)+ x(2:end,:,1,3);
					y(1,:,2:end) = y(1,:,2:end) - x(1,:,1:end-1,3) + x(1,:,2:end,3);
					y(1,:,1)=y(1,:,1)+x(1,:,1,3);
					% yy
					y(:,1,:)=y(:,1,:) + x(:,1,:,4);
					y(:,2,:)=y(:,2,:) + x(:,2,:,4) - 2*x(:,1,:,4);
					y(:,3:end,:)=y(:,3:end,:) + x(:,1:end-2,:,4) - 2*x(:,2:end-1,:,4) + x(:,3:end,:,4);
					% yz
					y(:,2:end,2:end)=y(:,2:end,2:end) + x(:,1:end-1,1:end-1,5) - x(:,1:end-1,2:end,5) - x(:,2:end,1:end-1,5) + x(:,2:end,2:end,5);
					y(:,2:end,1) = y(:,2:end,1) -x(:,1:end-1,1,5)+ x(:,2:end,1,5);
					y(:,1,2:end) = y(:,1,2:end) - x(:,1,1:end-1,5) + x(:,1,2:end,5);
					y(:,1,1)=y(:,1,1)+x(:,1,1,5);
					% zz
					y(:,:,1)=y(:,:,1) + x(:,:,1,6);
					y(:,:,2)=y(:,:,2) + x(:,:,2,6) - 2*x(:,:,1,6);
					y(:,:,3:end)=y(:,:,3:end) + x(:,:,1:end-2,6) - 2*x(:,:,2:end-1,6) + x(:,:,3:end,6);
					
        	end
		end
		 function y = HtH(this,x) %  Apply the HtH matrix
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizein);
            nidx = 0;
            y = zeros(this.sizein);
            % switch according to the number of dimension of the input
            switch(this.ndms)
                % 2 dimension
                case(2)
                	% xx
					y(1,:)= y(1,:) + x(3,:) -2*x(2,:) + x(1,:);
					y(2,:)=y(2,:) + x(4,:) - 4*x(3,:) + 5*x(2,:) - 2*x(1,:);
					y(3:end-2,:)=y(3:end-2,:) + 6*x(3:end-2,:) - 4*x(4:end-1,:) + x(5:end,:) - 4*x(2:end-3,:) + x(1:end-4,:);
					y(end-1,:)=y(end-1,:) + 6*x(end-1,:) - 4*x(end-2,:) - 4*x(end,:) + x(end-3,:);
					y(end,:)=y(end,:)+6*x(end,:)-4*x(end-1,:) + x(end-2,:);
					% xy
					y(1,1)=y(1,1) + x(2,2) - x(2,1) - x(1,2) + x(1,1);
					y(2:end-1,2:end-1)=y(2:end-1,2:end-1) + 4*x(2:end-1,2:end-1) - 2*x(2:end-1,1:end-2) - 2*x(1:end-2,2:end-1) -2*x(2:end-1,3:end) - 2*x(3:end,2:end-1) + x(1:end-2,1:end-2) + x(1:end-2,3:end) + x(3:end,1:end-2) + x(3:end,3:end); 
					y(1,2:end-1) = y(1,2:end-1) - 2*x(2,2:end-1) + 2*x(1,2:end-1) - x(1,3:end) + x(2,1:end-2) + x(2,3:end) - x(1,1:end-2);
					y(2:end-1,1) = y(2:end-1,1) - 2*x(2:end-1,2) + 2*x(2:end-1,1) - x(3:end,1) + x(1:end-2,2) + x(3:end,2) - x(1:end-2,1);
					y(end,2:end-1)=y(end,2:end-1) + 4*x(end,2:end-1) - 2*x(end,1:end-2) - 2*x(end,3:end) -2*x(end-1,2:end-1) + x(end-1,1:end-2) + x(end-1,3:end);
					y(2:end-1,end)=y(2:end-1,end) + 4*x(2:end-1,end) - 2*x(1:end-2,end) - 2*x(3:end,end) -2*x(2:end-1,end-1) + x(1:end-2,end-1) + x(3:end,end-1);
					y(1,end)=y(1,end) + 2*x(1,end) + x(2,end-1) -2*x(2,end) - x(1,end-1);
					y(end,1)=y(end,1) + 2*x(end,1) + x(end-1,2) -2*x(end,2) - x(end-1,1);
					y(end,end)=y(end,end) + 4*x(end,end) - 2*x(end,end-1) -2*x(end-1,end) + x(end-1,end-1);
                	% yy
					y(:,1)= y(:,1) + x(:,3) -2*x(:,2) + x(:,1);
					y(:,2)=y(:,2) + x(:,4) - 4*x(:,3) + 5*x(:,2) - 2*x(:,1);
					y(:,3:end-2)=y(:,3:end-2) + 6*x(:,3:end-2) - 4*x(:,4:end-1) + x(:,5:end) - 4*x(:,2:end-3) + x(:,1:end-4);
					y(:,end-1)=y(:,end-1) + 6*x(:,end-1) - 4*x(:,end-2) - 4*x(:,end) + x(:,end-3);
					y(:,end)=y(:,end)+6*x(:,end)-4*x(:,end-1) + x(:,end-2);
                % 3 dimensions
                case(3)
            end
    	end
    end
end

