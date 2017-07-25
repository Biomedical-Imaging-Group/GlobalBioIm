classdef (Abstract) Map < handle
  % comments
  %
  
	properties
	name = 'none'           % name of the linear operator

	isinvertible = false;   % true if the operator is invertible
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
	  if ~this.memoizeOpts.apply || ~isequal(this.memoCache.apply.in, x)
		y = this.apply_(x);
		if this.memoizeOpts.apply
		  this.memoCache.apply.in = x;
		  this.memoCache.apply.out = y;
		end
	  else
		y = this.memoCache.apply.out;
	  end
		
	  % check output size
	  if ~Map.checkSize(y, this.sizeout)
		warning('Output of apply was size [%s], didn''t match stated sizeout: [%s].',...
		  num2str(size(y)), num2str(this.sizeout));
	  end
	end
	
	function applyJacob(this, x, v)
	  
	end
	
	
	
  end
  
  methods (Access=protected)
        function y = apply_(this, x) 
        	% **(Abstract method)** Apply the linear operator
        	%
        	% :param x: \\(\\in X\\)
        	% :returns y: \\(= \\mathrm{Hx}\\)
        	
            error('Apply method not implemented');
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