classdef (Abstract) Map < handle
  % comments
  %
  
  properties
	  name = 'none'           % name of the linear operator
	  
	  isInvertible = false;   % true if H.applyInverse(  ) will work %todo fix capitalization everywhere
	  isDifferentiable = false; % true if H.applyJacobianT(   ) will work 
	  
	  % todo: decide on how to handle checking for inverse, grad, ... being
	  % implemented and working (no divide by zero).
	  %implementedMetods = struct('applyJacobianT', false);
	  
	  isComplex = false;      % true is the operator is complex %todo fix capitalization everywhere
	  
	  
	  norm=-1;                % norm of the operator
	  
	  sizein;                 % dimension of the right hand side vector space
	  sizeout;                % dimension of the left hand side vector space
	  
	  memoizeOpts = struct('apply', false);
	  doPrecomputation = false;
	  
	
	end
	
	properties (SetAccess = private)
	  memoCache = struct('apply', struct('in', [], 'out', []));
	  precomputeCache = struct();
	end

  methods (Sealed)
	
	function y = apply(this, x)
	  % check input size
	  if ~checkSize(x, this.sizein) 
		error('Input to apply was size [%s], didn''t match stated sizein: [%s].',...
		  num2str(size(x)), num2str(this.sizein));
	  end
	  
	  % memoize
	  y = this.memoize('apply', @this.apply_, x);
		
	  % check output size
	  if ~checkSize(y, this.sizeout)
		warning('Output of apply was size [%s], didn''t match stated sizeout: [%s].',...
		  num2str(size(y)), num2str(this.sizeout));
	  end
	end
	
	function x = applyJacobianT(this, y, v)
	  % x = J(v)^*{y}
	  
	  % check input size
	  if ~checkSize(y, this.sizeout)
		error('Input of y applyJacobianT was size [%s], didn''t match stated sizeout: [%s].',...
		  num2str(size(y)), num2str(this.sizeout));
	  end
	  
	  if ~checkSize(v, this.sizein)
		error('Input to v applyJacobianT was size [%s], didn''t match stated sizein: [%s].',...
		  num2str(size(v)), num2str(this.sizein));
	  end
	  
	  % memoize
	  x = this.memoize('applyJacobianT', @this.applyJacobianT_, {y, v});

	  % check output size
	  if ~checkSize(x, this.sizein)
		warning('Output of applyJacobianT was size [%s], didn''t match stated sizein: [%s].',...
		  num2str(size(x)), num2str(this.sizein));
	  end
	  
	  
	end
	
	function G = makeComposition(this, H)
	  % TODO: check conformable size
	  G = this.makeComposition_(this, H);
	  
	end
	
	function x = inverse(this, y)
		% TODO: size checking
		
		% **(Abstract method)** Apply \\(\\mathrm{H}^{-1}\\) (if applicable)
		%
		% :param x: \\(\\in Y\\)
		% :returns y: \\(= \\mathrm{H^{-1}x}\\)
		%
		
		if this.isinvertible
			error('inverse not implemented');
		else
			error('Operator not invertible');
		end
	end


	
  end
  
  
  methods (Access = protected)
	  function y = apply_(this, x)
		  error('apply method not implemented');
	  end
	  
	  function x = applyJacobianT_(this, y, v)
		  error('Jacobian method not implemented');
	  end
	  
	  
	  function x = applyInverse_(this, y)
		  error('Inverse method not implemented');
	  end
	  
	  % composition methods
	  function G = makeComposition_(this, H)
		  G = MapComposition({this, H});
	  end
	  

	  
		
		
		% utility methods		
		function y = memoize(this, fieldName, fcn, xs)
		  
		  if ~iscell(xs) % handle single input case
			xs = {xs};
		  end
		  
		  if ~this.memoizeOpts.(fieldName) || ~isequal(this.memoCache.(fieldName).in, xs)
			y = fcn(xs{:});
			if this.memoizeOpts.(fieldName)
			  this.memoCache.(fieldName).in = xs;
			  this.memoCache.(fieldName).out = y;
			end
		  else
			y = this.memoCache.(fieldName).out;
		  end
		end
		
		
  end

end