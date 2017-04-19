classdef MulLinOp < LinOp
    %% MulLinop : Multiplication of linear operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = SumLinop(LinOp1,LinOp2)
    % Multiplication of LinOps:
    % Obj = LinOp1 * LinOp2
    %
    %
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
        LinOp1
        LinOp2
        isnum
		isHTH = 0;
		isHHt = 0;
    end
    
    methods
        function this = MulLinOp(LinOp1, LinOp2)
            
            this.name ='MulLinOp';
            this.LinOp1 = LinOp1;
            this.LinOp2 = LinOp2;
            
			% strcmp is different than isa because it doesn't check all
			% superclasses as well
			if strcmp(class(LinOp1), 'adjoint') && LinOp1.TLinOp == LinOp2
				this.isHTH = true;
				this.issquare = true;
			elseif strcmp(class(LinOp2), 'adjoint') && LinOp2.TLinOp == LinOp1
				this.isHHt = true;
				this.issquare = true;
			end

			if isnumeric(LinOp1) && LinOp2.norm ~= -1
				this.norm = LinOp1 * LinOp2.norm;
			elseif ~isnumeric(LinOp1) && LinOp1.norm ~= -1 && this.LinOp2.norm ~= -1
				this.norm = this.LinOp1.norm * this.LinOp2.norm;
            else
                this.norm=-1;
			end
					
		
            if isnumeric(LinOp1)
                this.isnum =1;
                if (~isreal(LinOp1)) || LinOp2.iscomplex
                    this.iscomplex= true;
                else
                    this.iscomplex= false;
                end
                
                if all(LinOp1) && LinOp2.isinvertible
                    this.isinvertible= true;
                else
                    this.isinvertible= false;
                end
                
                this.issquare= LinOp2.issquare;             
                
                %                 if isscalar(LinOp1)
                %              LinOp1 = Scaling(LinOp1);
                %                 else
                %              LinOp1 = Diagonal(LinOp1);
                %                 end
            else
                assert(isa(LinOp1,'LinOp'),'MulLinOp: First input should be a LinOp');
                assert(isa(LinOp2,'LinOp'),'MulLinOp: Second input should be a LinOp');
                
                assert(isempty(LinOp1.sizein) || isequal(LinOp1.sizein, LinOp2.sizeout),'size of LinOp not conformable');
                this.sizein = LinOp2.sizein;
                this.sizeout = LinOp1.sizeout;
                if LinOp1.iscomplex || LinOp2.iscomplex
                    this.iscomplex= true;
                else
                    this.iscomplex= false;
                end
                
                if LinOp1.isinvertible && LinOp2.isinvertible
                    this.isinvertible= true;
                else
                    this.isinvertible= false;
                end
                
                if LinOp1.issquare && LinOp2.issquare
                    this.issquare= true;
                else
                    this.issquare= false;
                end               
            end                 
        end
        
		function y = apply(this,x) % apply the operator
			if this.isHTH
				y = this.LinOp2.HtH(x);
			elseif this.isHHt
				y = this.LinOp1.HHt(x);
			elseif this.isnum
				y = this.LinOp1.*this.LinOp2.apply(x);
			else
				y = this.LinOp1.apply( this.LinOp2.apply(x));
			end
		end
		function y = adjoint(this,x) % apply the adjoint
			if this.isHTH || this.isHHt
				y = this.apply(x); % because self-adjoint
			elseif this.isnum
				y = this.LinOp2.adjoint(this.LinOp1.*x);
			else
				y = this.LinOp2.adjoint(this.LinOp1.adjoint(x));
			end
		end
		function y = HtH(this,x)
			if this.isHTH || this.isHTH 
				y = this.apply(this.apply(x)); % because self-adjoint
			elseif this.isnum
				y = this.LinOp2.adjoint(this.LinOp1.^2.*( this.LinOp2.apply(x)));
			else
				y = this.LinOp2.adjoint(this.LinOp1.HtH( this.LinOp2.apply(x)));
			end
        end
		function y = HHt(this,x)
			if this.isHTH || this.isHTH 
				y = this.apply(this.apply(x)); % because self-adjoint
			elseif this.isnum
				y = this.LinOp1.*(this.LinOp2.HHt( this.LinOp1.*x));
			else
				y = this.LinOp1.apply(this.LinOp2.HHt( this.LinOp1.adjoint(x)));
			end
		end
		function y = inverse(this,x) % apply the inverse
			if this.isinvertible
				y = this.LinOp2.inverse(this.LinOp1.inverse(x));
			else
				error('Operator not invertible');
			end
		end
        function y = adjointInverse(this,x) % apply the inverse
            if this.isinvertible             
                y = this.LinOp2.adjointInverse(this.LinOp1.adjointInverse(x));
            else
                error('Operator not invertible');
            end
        end
    end
end

