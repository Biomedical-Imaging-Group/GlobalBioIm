classdef LinOpXRay < LinOp
   % LinOpXray: 2D, parallel x-ray operator based on MATLAB's ``radon`` function.
   % 
   %
   % :param x: structure defining the reconstruction domain (always 2D)
   % :param x.size: 1x2 vector giving dimensions of the domain in pixels, must be of the form ``N*ones(1,2)``.
   % :param x.step: 1x2 vector giving the spacing of the pixels
   %                (in an apropirate unit of length, e.g. milimeteres).
   % :param thetas: 1x(number of projections) vector giving the projection directions in radians.
%                   Specifically, ``[cos(theta) sin(theta)]`` points in the direction of the positive axis in the projection coordinate system.
   %
   % **Note** We take rows as $$x$$ and columns as $$y$$ in the reconstruction domain.
   % **Note** The center of the reconstruction (and rotation axis) is assumed to be at ``floor((x.size+1)/2)``. 
   %
   % **Note** The center of each projection is at ``ceil(sizeout/2)``.
   %
   % **Note** Be warned that the adjoint is only accurate to around 20 dB for
   % for images of size 500x500 px. The error is smaller for larger
   % images.

   %% Copyright (C) 2019
   %  M. McCann michael.thompson.mccann@gmail.com
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
			g = iradon(c, this.thetas/pi*180 - 90, 'linear', 'none', 1, this.sizein(1));
			g = g * this.x.step(1);	
		end

		function y = applyHtH_(this,x) % apply the HtH matrix
			y = this.applyAdjoint(this.apply(x)); 
		end

	end


end

