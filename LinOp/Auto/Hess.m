function obj = Hess(varargin) 
warning('Hess is deprecated, please use LinOpHess'); 
obj = LinOpHess(varargin{:}); 
