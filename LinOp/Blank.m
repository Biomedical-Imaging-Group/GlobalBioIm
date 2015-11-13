classdef Blank < LinOp
	properties (SetAccess = public, GetAccess = public)
		
		
	end
	
	methods
		function this = Blank()
			
		end
		
		function Apply(this,~) % Apply the operator
			error('Operator not implemented');
		end
		function Adjoint(this,~) % Apply the adjoint
			error('Adjoint not implemented');
		end
		function y = HtH(this,x) %  Apply the HtH matrix
			y = this.Adjoint(this.Apply(x));
		end
		function y = HHt(this,x) %  Apply the HHt matrix
			if this.issquare   % HtH =HHt
				y = this.HtH(x);
			else
				y = this.Apply(this.Ajoint(x));
			end
		end
		function y = HtWH(this,x,W) %  Apply the HtH matrix
			if (isscalar(W) && isreal(W))
				y = this.HtH(x);
			else
				assert(isa(W,'LinOp'),'W must be a LinOp');
				y = this.Adjoint(W.Apply(this.Apply(x)));
			end
		end
		function Inverse(this,~) % Apply the inverse
			if this.isinvertible
				error('Inverse not implemented');
			else
				error('Operator not invertible');
			end
		end
		function AdjointInverse(this,~) % Apply the inverse
			if this.isinvertible
				error('AdjointInverse not implemented');
			else
				error('Operator not invertible');
			end
		end
	end
	
	methods (Private)
			
	end
end