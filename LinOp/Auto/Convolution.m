function obj = Convolution(varargin) 
warning('Convolution is deprecated, please use LinOpConv'); 
obj = LinOpConv(varargin{:}); 
