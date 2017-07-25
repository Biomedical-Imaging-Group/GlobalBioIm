classdef (Abstract) Map < handle
  % comments
  %
  
	properties
	name = 'none'           % name of the linear operator

	isinvertible = false;   % true if the operator is invertible
	
	% todo: decide on how to handle checking for inverse, grad, ... being
	% implemented and working (no divide by zero).
	%implementedMetods = struct('applyJacobianT', false);
	
	iscomplex = false;      % true is the operator is complex
	
	
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
	  if ~Map.checkSize(x, this.sizein) 
		error('Input to apply was size [%s], didn''t match stated sizein: [%s].',...
		  num2str(size(x)), num2str(this.sizein));
	  end
	  
	  % memoize
	  y = this.memoize('apply', @this.apply_, x);
		
	  % check output size
	  if ~Map.checkSize(y, this.sizeout)
		warning('Output of apply was size [%s], didn''t match stated sizeout: [%s].',...
		  num2str(size(y)), num2str(this.sizeout));
	  end
	end
	
	function x = applyJacobianT(this, y, v)
	  % x = J(v)^*{y}
	  
	  % check input size
	  if ~Map.checkSize(y, this.sizeout)
		error('Input of y applyJacobianT was size [%s], didn''t match stated sizeout: [%s].',...
		  num2str(size(y)), num2str(this.sizeout));
	  end
	  
	  if ~Map.checkSize(v, this.sizein)
		error('Input to v applyJacobianT was size [%s], didn''t match stated sizein: [%s].',...
		  num2str(size(v)), num2str(this.sizein));
	  end
	  
	  % memoize
	  x = this.memoize('applyJacobianT', @this.applyJacobianT_, {y, v});

	  % check output size
	  if ~Map.checkSize(x, this.sizein)
		warning('Output of applyJacobianT was size [%s], didn''t match stated sizein: [%s].',...
		  num2str(size(x)), num2str(this.sizein));
	  end
	  
	  
	end
	
	function G = compose(this, H)
	  % TODO: check conformable size
	  
	  G = this.compose_(H);
	  
	end
	
	  function y = inverse(this, x) 
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
  
  methods (Access = private)
        function y = apply_(this, x) 
        	% **(Abstract method)** Apply the linear operator
        	%
        	% :param x: \\(\\in X\\)
        	% :returns y: \\(= \\mathrm{Hx}\\)
        	
            error('Apply method not implemented');
		end
		
		function x = applyJacobianT_(this, y, v)
		     error('Jacobian method not implemented');
		end
		
		function G = this.compose_(this, H)
		  % todo: G = MapComposition({this, H});
		end
		
		function x = this.inverse_(this, y)
		  error('inverse method not implemented');
		end
  end
  
  methods (Access = protected)
		
		
		% utility functions		
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
  
 	methods (Static)
	  function sizeOK = checkSize(x, sz)
		  xSize = size(x);
		  if length(xSize) < length(sz)
			xSize(end+1:length(sz)) = 1;
		  elseif length(sz) < length(xSize)
			sz(end+1:length(xSize)) = 1;
		  end
		  
		  sizeOK = isequal(xSize,sz);
		  
		end
	end
end