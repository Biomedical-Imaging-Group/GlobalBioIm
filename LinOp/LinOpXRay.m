classdef LinOpXRay < LinOp
   % LinOpXray: 2D, parallel x-ray operator based on MATLAB's radon.
   % Be warned that the adjoint is only accurate to around 20 dB for
   % for images of size 500x500 px. The error is smaller for larger
   % images.
   %
   % :param x: structure defining the reconstruction domain (always 2D)
   %        x.size - 1x2 vector giving dimensions of the domain in pixels,
   %                 must be of the form N*ones(1,2)
   %        x.step - 1x2 vector giving the spacing of the pixels (in
   %                 an apropirate unit of length, e.g. milimeteres)
   %        the center is assumed to be at floor((x.size+1)/2)
   %
   %  the center of each projection is at ceil(sizeout/2)
 
    properties
     x
	 y
	 thetas 
    end
    
    methods

		function this = LinOpXRay(x, thetas)
            this.name ='XRay';
			
			% input checking for x
			if length(unique(x.size))>1
				error('only square x domains are supported');
			end
			
			if length(x.size) ~= 2
				error('only 2d supported');
			end
			
			if length(unique(x.step))>1
				error('only isotropic x domains are supported');
			end
			
			this.x = x;
			this.sizein = x.size;

			% setting output size
			y.size = 2*ceil(norm(x.size-floor((x.size-1)/2)-1))+3;
            this.sizeout= [y.size length(thetas)];

			this.y = y;
			this.thetas = thetas;
			
			% make copies of private image processing toolbox functions so
			% that we can call them directly
			if isempty(which('iradonmex'))
				filename = which('iradonmex', '-all');
				privateDir = fullfile(fileparts(mfilename('fullpath')), 'private');
				if ~exist(privateDir, 'dir')
					mkdir(privateDir)
				end
				copyfile(filename{1}, privateDir)
			end

			% estimate norm, assuming shift invariance
			d = zeros(x.size);
			d(ceil(end/2), ceil(end/2)) = 1;
			
			HTHd = fft2(this.applyAdjoint_(this.apply_(d)));

			this.norm = sqrt(max(abs(HTHd(:))));

			
		end
		
	
	end
	
	methods (Access = protected)
		
		function y=apply_(this,x) % apply the operator
			y=radon(x,this.thetas/pi*180 - 90, this.y.size) * this.x.step(1);
		end

		function g = applyAdjoint_(this,c) % apply the adjoint
			xVects = (1:this.x.size(1)) - floor((this.x.size(1)+1)/2);
			[X0, X1] = ndgrid(xVects);
			g = iradonmex(size(X0,1), this.thetas - pi/2, X1./this.x.step(2), -X0./this.x.step(1), c, 1); % 1 - linear interp
			g = g * this.x.step(1);	
		end

		function y = applyHtH_(this,x) % apply the HtH matrix
			y = this.applyAdjoint(this.apply(x)); 
		end

	end


end

