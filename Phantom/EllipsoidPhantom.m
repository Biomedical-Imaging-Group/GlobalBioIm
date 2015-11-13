classdef EllipsoidPhantom < XRayPhantom
  properties
	E; 
	D; % number of dimensions
  end
  
  methods(Static)
	function f = example(num)
	  % example of the usage of the ellipsoid phantom
	  if ~exist('num', 'var') || isempty(num)
		num = 1;
	  end
	  
	  switch num
		%----------------------------------------------------------------
		case 1
		  % each row is an ellipse.
		  % intensity, a, b, centerX, centerY, rotation, just like 'help
		  % phantom'
		  ells = [1, 10, 3, 11, 5, pi/8
			.4, 5, 5, 0, 0, 0
			2, 10, 7, -1, -10, 2*pi/3
			2 .5 .5 -5 7 0]
		  p = EllipsoidPhantom(ells);
		  
		  xStart = [-25 -25];
		  xStep = [.25 .25];
		  xSize = [200 200];
		  f = p.eval(xStart, xStep, xSize);
		  
		  figure
		  % rows of f correspond to x_0, cols to x_1 , pages to x_2, so we need to
		  % transpose to plot it in the conventional way, with "x" going
		  % left-right and "y" going up-down
			imagesc( (0:xSize(1)-1) * xStep(1) + xStart(1), ...
				(0:xSize(2)-1) * xStep(2) + xStart(2), ...
				f');
			axis xy
			colorbar
			xlabel('x_0')
		  ylabel('x_1')
		%----------------------------------------------------------------  
		case 2 % x-ray projections of a simple 2D phantom
		  ells = [1, 10, 3, 11, 5, pi/10 + pi/2]
		  p = EllipsoidPhantom(ells);
		  
		  
		  %% calculate a sinogram
		  numProj = 10; % number of projection angles
		  P = zeros(1, 2, numProj);
		  t = (0:numProj-1) * pi/numProj;
		  for i = 1:numProj
				P(:,:,i) = [cos(t(i)) sin(t(i))];
			end
			
		  yStart = -25;
		  yStep = 1;
		  ySize = 50;
		  g = p.xray(yStart, yStep, ySize, P);
		  
		  sino_h = figure;
		  imagesc(g, 'XData', t, 'YData', (0:ySize-1)*yStep + yStart)
		  xlabel('perpendicular to projection angle')
		  ylabel('y')
		  axis manual;
		  title('sinogram');
		  
		  %% display the phanotm
		  xStart = [-25 -25];
		  xStep = [.25 .25];
		  xSize = [200 200];
		  f = p.eval(xStart, xStep, xSize);
		  
		  f_h = figure;
		  % rows of f correspond to x_0, cols to x_1 , pages to x_2, so we need to
		  % transpose to plot it in the conventional way, with "x" going
		  % left-right and "y" going up-down
		  imagesc((0:xSize(1)-1) * xStep(1) + xStart(1), ...
			(0:xSize(2)-1) * xStep(2) + xStart(2), ...
			f');
		  axis xy
		  colorbar
		  axis manual % so the axes don't change when we plot lines
		  xlabel('x_0')
		  ylabel('x_1')
		  
		  %% display some lines

		  styles = {'g', 'y', 'r'};

		  inds = [2 7 7
			38 18 37];
		  for i = 1:size(inds,2)
				pInd = inds(1,i);
				yInd = inds(2,i);
				
				figure(f_h)
				hold on
				y = calcGridPoints(yStart, yStep, ySize);
				x = P(:,:,pInd)' * y;
				u = null( P(:,:,pInd) );
				%plot(x(1,yInd), x(2,yInd), 'o')
				l = [x(:,yInd) + xStep(1)*xSize(1)*u, x(:,yInd) - xStep(2)*xSize(2)*u];
				plot(l(1,:), l(2,:), styles{i});
				
				figure(sino_h)
				hold on
				plot(t(pInd), y(yInd), ['.' styles{i}], 'markersize', 20);
			end
			
	
		%----------------------------------------------------------------  
		case 3 % 3-d example
		  % each row is an ellipse.
		  ells = [1, 5, 3, 9, 0, 0, 0, 0, 0, 0
			1, 7, 17, 3, 10, 3, 0, 0, 0, 0
			3, 10, 1, 1, -10, -5, 0, pi/4, 0, 0
			2, 10, 4, 20, -15, 10, 0, pi/8, 2*pi/3, 0]
		  
		  p = EllipsoidPhantom(ells);
		  
		  xStart = [-25 -25 -25];
		  xStep = [.25 .25 .25];
		  xSize = [200 200 200];
		  f = p.eval(xStart, xStep, xSize);
		  
		  
			for slice = [65, 101, 121];
				figure
				imagesc((0:xSize(1)-1) * xStep(1) + xStart(1), ...
					(0:xSize(2)-1) * xStep(2) + xStart(2), ...
					f(:,:,slice)', [0 5]);
				axis xy
				colorbar
				xlabel('x_0')
				ylabel('x_1')
				title(sprintf('x_2 = %g', xStep(3)*(slice-1) + xStart(3)));
			end
			
		  %% take x-ray measurements
		  % setup the y grid, which is 2D
		  yStart = [-30 -30];
		  yStep = [.25 .25];
		  ySize = [300 300];
		  
		  % set angles determining the projection directions
		  phi =   [0  0     0     pi/4  pi/2  pi/4];
		  theta = [0  -pi/2  pi/2  0     0     pi/4];
		  psi =   [0  0     0     0     0     0   ];
		  numProj = length(phi); % number of projection angles
		  
		  % generate a basis orthogonal to the projection direction for
		  % each projection.
		  [r0, r1, r2] = Euler3D(phi, theta, psi);
		  P = reshape([r0(:)'; r1(:)'], 2, 3, numProj);
		  
		  % take the measurements
		  g = p.xray(yStart, yStep, ySize, P);
		  
		  % plot
			
			
			for gInd = [1 2];
				figure
				imagesc( (0:ySize(1)-1)*yStep(1) + yStart(1), ...
					(0:ySize(2)-1)*yStep(2) + yStart(2), ...
					g(:,:,gInd)');
				axis xy
				colorbar
				xlabel(sprintf('y_0 = [%.0f %.0f %.0f]', r0(:,gInd)));
				ylabel(sprintf('y_1 = [%.0f %.0f %.0f]', r1(:,gInd)));
				title(sprintf('projection along [%.0f %.0f %.0f]', r2(:,gInd)));
			end
			
			figure
			volume_plot(xStart, xStep, xSize, f, .5);
			
		end % switch
	  
	  
	end % example func
	end % static methods
  
  methods %----------------------------------------------------------------
	
	function obj = EllipsoidPhantom(E)
		% obj = EllipsoidPhantom(E)
		%
		% E - defines the ellipsiods in the phantom.
		%     FOR A 2D PHANTOM: E is P X 6, where P is the number of ellipses
		%     E(p, 1) = intensity
		%     E(p, 2) = semi-axis 1 length
		%     E(p, 3) = semi-axis 2 length
		%     E(p, 4) = center coordinate 1
		%     E(p, 5) = center coordinate 2
		%     E(p, 6) = angle (radians) between axis 1 and the horizontal
		%
		%     FOR A 3D PHANTOM: E is P X 10
		%     E(p, 1) = intensity
		%     E(p, 2) = semi-axis 1 length
		%     E(p, 3) = semi-axis 2 length
		%     E(p, 4) = semi-axis 2 length
		%     E(p, 5) = center coordinate 1
		%     E(p, 6) = center coordinate 2
		%     E(p, 7) = center coordinate 2
		%     E(p, 8:10) = phi, theta, psi giving the z-y-z euler angles for
		%     the semi-axes
		
	  obj.E = E;
	  if size(E, 2) == 6
		obj.D = 2;
	  elseif size(E, 2) == 10
		obj.D = 3;
	  else
		error('unknown format for input E');
	  end
	end
	
	function f = eval(obj, xStart, xStep, xSize)
		% f = eval(obj, xStart, xStep, xSize)
		
	  if length(xStart) ~= obj.D || ...
		  length(xStep) ~= obj.D || ...
		  length(xSize) ~= obj.D
		error('dimension mismatch');
	  end
	  
	  % setup output grid
	  x = makeGrid(xStart, xStep, xSize);
	  f = zeros(xSize);
	  
	  % do calcs for each ellipse
	  for eInd = 1:size(obj.E,1);
		% find the rotation matrix, Phi'
		if obj.D == 2
		  R = [cos(obj.E(eInd, 6)) -sin(obj.E(eInd, 6))
			sin(obj.E(eInd, 6)) cos(obj.E(eInd, 6))]';
		else 
		  phi = obj.E(eInd, 1+2*obj.D+1);
		  theta = obj.E(eInd, 1+2*obj.D+2);
		  psi = obj.E(eInd, 1+2*obj.D+3);
		  [r1, r2, r3] = Euler3D(phi, theta, psi);
		  R = [r1 r2 r3]';
		end
		
		% handle the center offset
		offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
		pos = R * bsxfun(@minus, x, offset);
		
		% handle the rotation
		axes = obj.E(eInd, (1:obj.D) + 1)';
		vals = sum( bsxfun(@times, pos.^2, 1./axes.^2) ); 
		
		% add ellipse to the output
		f(vals <= 1) = f(vals <= 1) + obj.E(eInd,1);
	  end
	  
	  
	  
	end
	
	function g = xray(obj, yStart, yStep, ySize, P)
		% g = xray(obj, yStart, yStep, ySize, P)
		
	  if length(yStart) ~= obj.D-1 || ...
		  length(yStep) ~= obj.D-1 || ...
		  length(ySize) ~= obj.D-1
		error('dimension mismatch');
	  end
	  
	  
	  y = makeGrid(yStart, yStep, ySize);
	  g = zeros([ySize size(P, 3)]);
	  
	  % for each ellipse...
	  for eInd = 1:size(obj.E,1);
		% get the center and axis lengths
		offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
		axes = obj.E(eInd, (1:obj.D) + 1)';
		
		% find the rotation matrix for the ellipse
		if obj.D == 2
		  R = [cos(obj.E(eInd, 6)) -sin(obj.E(eInd, 6))
			sin(obj.E(eInd, 6)) cos(obj.E(eInd, 6))]';
		else
		  phi = obj.E(eInd, 1+2*obj.D+1);
		  theta = obj.E(eInd, 1+2*obj.D+2);
		  psi = obj.E(eInd, 1+2*obj.D+3);
		  [r1, r2, r3] = Euler3D(phi, theta, psi);
		  R = [r1 r2 r3]';
		end
		
		% for each projection angle
		for pInd = 1:size(P,3)
		  % find the intersection of the ellipse with the line 
		  % y(t) = x + t*u
		  
		  % we do this in the coordinate system of the ellipse, so that it
		  % is axis-aligned.
		  x = R * bsxfun(@minus, P(:,:,pInd)' * y, offset);
		  u = R * null( P(:,:,pInd) );
		  
		  % find real roots of the quadratic defined by the intersection
		  % Ax^2 + Bx + C = 0
		  A = sum(u.^2 ./ axes.^2);
		  B = 2 * sum( bsxfun(@times, u./axes.^2, x) );
		  C = sum( bsxfun(@times, x.^2, 1 ./ axes.^2) ) - 1;
		  
		  d = sqrt( B.^2 - 4 * A * C );
		  r1 = ( -B + d )/ (2*A);
		  r2 = ( -B - d )/ (2*A);
		  
		  h = imag(r1) == 0 & imag(r2) == 0;
		  
		  hInd = find(h)' + prod(ySize) * (pInd-1)';
		  g(hInd) = g(hInd) + abs(r1(h) - r2(h))' * obj.E(eInd,1);
		end
	  end
	  
	end
  end
  
  
end