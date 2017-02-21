function obj = Grad(varargin) 
warning('Grad is deprecated, please use LinOpGrad'); 
obj = LinOpGrad(varargin{:}); 
