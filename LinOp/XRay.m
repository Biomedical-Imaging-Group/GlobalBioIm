classdef XRay < LinOp
	properties (SetAccess = public, GetAccess = public)
		D
		
		xStart
		xStep
		xSize
		
		yStart
		yStep
		ySize
		Ps
		
		phi
		numProj
		
		yUpsampleRate = 3;
		HTHHat
		
		
		
	end
	
	methods
		function this = XRay(xStart, xStep, xSize, yStart, yStep, ySize, Ps, phi)
			this.D = size(xStart,2);
			this.phi = phi;
			this.Ps = Ps;
			this.numProj = size(Ps,3);
			
			
			this.xStart = xStart;
			this.xStep = xStep;
			this.xSize = xSize;
			this.sizein = xSize;
					
			
			this.yStart = yStart;
			this.yStep = yStep;
			this.ySize = ySize;
			this.sizeout = [ySize this.numProj];
		end
		
		function Apply(this,~) % Apply the operator
			error('Operator not implemented');
		end
		
		function HTg = Adjoint(this,g) % Apply the adjoint
			HTg = this.adjoint_LUT(g, this.yStart, this.yStep, this.ySize, this.Ps, ...
				this.xStart, this.xStep, this.xSize, this.phi, this.yUpsampleRate);
		end
		
		function g = HtH(this, c)
			if isempty(this.HTHHat)
				phiRadiusPixels = this.phi.spaceRadius / this.xStep(1);
				
				kSize(1:2) = this.xSize(1:2)*2 - 1;
				if this.D==3
					kSize(3) = this.xSize(3) + (phiRadiusPixels*2+1) - 1;
				end
				kStep = this.xStep;
				kStart = this.xStart - ((kSize-this.xSize )/2 .* kStep);
				
				yUpStep =  ones(1, this.D-1)/2; % 2x upsampling
				doMex = this.D==3;
				isIsotropic = false;
				bigMode = false;
				
				h = this.HTH_kernel_LUT(	yUpStep, kStart, kStep, kSize, doMex, isIsotropic, bigMode);
				
				centerPos = (-kStart ./ kStep) + 1;
				h = circshift(h, -centerPos+1);
				
				this.HTHHat = fftn( h );
				clear h
				
			end
			g = applyKernel(this.HTHHat, c);
			
			
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
	
	methods (Access = private)
		function HTg = adjoint_LUT(this, g, yStart, yStep, ySize, Ps, ...
				xStart, xStep, xSize, phi, yUpsampleRate)
			% function to compute the x-ray adjoint cTilde = P'g in the space domain in
			% 2D or 3D via LUT
			%
			% input:
			% g - measurements
			%
			% Mike McCann 2015
			
			%% handle input
			D = size(Ps,2);
			if ~( D == 2 || D == 3)
				error('Ps suggests this is a %d-D problem, which we do not handle', D);
			end
			
			if isempty(phi.spaceRadius) || ~isfinite(phi.spaceRadius)
				error('LUT adjoint only works for phi''s that have finite support in space');
			end
			
			if length(unique(phi.isotropicGroups)) > 1
				error('LUT adjoint currently only handles fully isotropic phi''s');
			end
			
			numProj = size(Ps,3);
			
			%% setup for main loop
			% upsampling setup
			yUpStep = yStep ./ yUpsampleRate;
			yUpSize = ySize .* yUpsampleRate;
			
			upSampleVects = make_grid_vectors(ones(1,D-1), yUpsampleRate, ySize);
			yInds = make_grid_vectors(ones(1,D-1), ones(1,D-1), ySize);
			
			% convolution setup
			r = phi.spaceRadius ./yUpStep;
			phiGrid = makeGrid(-yUpStep .* ceil(r), yUpStep, 2*ceil(r)+1);
			phiVals = phi.xray(phiGrid);
			
			fftSize = yUpSize + 2*ceil(r)+1 - 1;
			phiHat = fftn( rot90(reshape(phiVals, [2*ceil(r)+1 1]), 2), [fftSize 1] );
			% rot90 handles reversing the mask, which I'm not sure is what we should be
			% doing, but it's consistent with previous code
			
			% interp setup
			%HTg = zeros(xSize);
			HTg = zeros(prod(xSize), 1);
			xkGrid = makeGrid(xStart, xStep, xSize)'; % might be big! :(
			
			yUpStart = yStart - ceil(r) .* yUpStep;
			gConvPhiVects = make_grid_vectors(yUpStart, yUpStep, fftSize);
			
			
			%% main loop
			gUp = zeros( yUpSize );
			% (we chose to process projInd by projInd instead of doing the upsamping and
			% convolution in a single step to reduce memory consumption)
			for projInd = 1:numProj
				% upsample g
				gUp( upSampleVects{:} ) = g( yInds{:}, projInd );
				
				% compute g conv phi
				gConvPhi = ifftn( fftn( gUp, [fftSize 1] ) .* phiHat, 'symmetric');
				
				% project the xGrid
				Px = xkGrid * Ps(:, :, projInd)';
				
				% interpolate
				F = griddedInterpolant(gConvPhiVects, gConvPhi);
				HTg = HTg + F(Px);
			end
			HTg = reshape(HTg, xSize);
		end
		
		
		function HTHkernel = HTH_kernel_LUT(this, yStep, xStart, xStep, xSize, doMex, isIsotropic, bigMode)
			% todo: this function is still a mess. in short, it should intelligently
			% make use of the isotropy and separability of the kernel
			%
			% input:
			%
			% Mike McCann 2015
			
			
			replaceVal = 0;
			numProj = this.numProj;
			D = this.D;
			phi = this.phi;
			Ps = this.Ps;
			
			if length(yStep) ~= this.D-1
				error('bad yStep');
			end
			
			
			%% compute g conv phi via filtering
			if ~isnan(phi.spaceRadius) % space limited
				r = phi.spaceRadius ./ yStep;
				phiGrid = makeGrid(-yStep .* ceil(r), yStep, 2*ceil(r)+1);
				phiVals = phi.xray(phiGrid);
				
				if D == 2
					rLUT = squeeze(imfilter(permute(phiVals, [2 3 1]), permute(phiVals, [2 3 1]), 'full', 'corr')) * prod(yStep);
				elseif D == 3
					rLUT = imfilter(reshape(phiVals, 2*ceil(r)+1), reshape(phiVals, 2*ceil(r)+1), 'full', 'corr') * prod(yStep);
				end
				
			end
			
			
			%% do the interpolation
			if doMex % fast, memory efficient MEX version (no ndgrid), but right now does not
				% make use of isotropy ----------------------------------------
				if D ~= 3; error('no 2d MEX version'); end
				yUpStart = -floor(size(rLUT)/2) .* yStep;
				HTHkernel = HTH_LUT_inner(rLUT, yUpStart, yStep, size(rLUT), Ps, xStart, xStep, xSize);
				
			elseif bigMode %---------------------------------------------------
				if D ~= 3; error('no 2d bigmode'); end
				
				HTHkernel = zeros(xSize);
				
				xkGrid = makeGrid(xStart(1:2), xStep(1:2), xSize(1:2));
				rLUT = rLUT(  (0:size(rLUT,1)-1) - (2*ceil(r(1))) == 0, :);
				points = ( (0:size(rLUT,2)-1) - (2*ceil(r(2))) ) * yStep(2);
				for xLevel = 1:xSize(3)
					fprintf('%d\n', xLevel);
					xkGrid(3,:) = xStart(3) + (xLevel-1)*xStep(3);
					
					if isIsotropic
						
						for thetaInd = 1:numProj
							P = Ps(:,:,thetaInd);
							offset = P * xkGrid;
							sOffset = sqrt( sum(offset.^2) );
							
							vals = interp1(points, rLUT(:), sOffset, 'pchip', replaceVal);
							
							HTHkernel(:, :, xLevel) = HTHkernel(:, :, xLevel) + reshape(vals, xSize(1), xSize(2));
						end
						
					else % not isotropic
						
						points = cell(D-1,1);
						for dInd = 1:D-1
							points{dInd} = ( (0:size(rLUT,dInd)-1) - (2*ceil(r(dInd))) ) * yStep(dInd);
						end
						
						
						
						for thetaInd = 1:numProj
							P = Ps(:,:,thetaInd);
							offset = P * xkGrid;
							
							% need to swap points 2 and 1 and offsets(2,:) and offset(1,:) because
							% of the rows = y problem
							vals = interp2(points{2}, points{1}, rLUT(:,:), offset(2,:), offset(1,:), 'linear', replaceVal); %'cubic'
							
							HTHkernel(:) = HTHkernel(:) + vals';
						end
					end
				end % xLevel loop
				
			else % not bigmode or mex ---------------------------------
				
				HTHkernel = zeros(xSize);
				xkGrid = makeGrid(xStart, xStep, xSize);
				
				if isIsotropic
					
					rLUT = rLUT(  (0:size(rLUT,1)-1) - (2*ceil(r(1))) == 0, :);
					points = ( (0:size(rLUT,2)-1) - (2*ceil(r(2))) ) * yStep(2);
					
					for thetaInd = 1:numProj
						P = Ps(:,:,thetaInd);
						offset = P * xkGrid;
						sOffset = sqrt( sum(offset.^2) );
						
						vals = interp1(points, rLUT(:), sOffset, 'pchip', replaceVal);
						
						HTHkernel(:) = HTHkernel(:) + vals';
					end
					
					
				else % not isotropic
					
					points = cell(D-1,1);
					for dInd = 1:D-1
						points{dInd} = ( (0:size(rLUT,dInd)-1) - (2*ceil(r(dInd))) ) * yStep(dInd);
					end
					
					
					for thetaInd = 1:numProj
						P = Ps(:,:,thetaInd);
						offset = P * xkGrid;
						
						if D == 2
							vals = interp1(points{1}, rLUT(:), offset, 'pchip', replaceVal);
						elseif D == 3
							% need to swap points 2 and 1 and offsets(2,:) and offset(1,:) because
							% of the rows = y problem
							vals = interp2(points{2}, points{1}, rLUT(:,:), offset(2,:), offset(1,:), 'linear', replaceVal); %'cubic'
						end
						HTHkernel(:) = HTHkernel(:) + vals';
					end
					
				end
			end % mex, bigmode, normal mode
		end
		
	end
	
	
	
end
