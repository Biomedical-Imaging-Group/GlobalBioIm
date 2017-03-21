function obj = Identity(varargin) 
warning('Identity is deprecated, please use LinOpIdentity'); 
obj = LinOpIdentity(varargin{:}); 
