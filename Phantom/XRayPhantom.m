classdef XRayPhantom
  % superclass for xray phantoms
  
  properties

	
  end
  
  methods(Abstract)
	g = xray(obj, yStart, yStep, ySize, P)
	f = eval(obj, xStart, xStep, xSize)
  end
  
  methods
	function g = measure(obj, edges, thetas)
		error('function is not written!')
	  % unfinished function to take xray measurements of an arbitrary
	  % phantom numerically
	  
	  for thetaInd = 1:length(thetas)
		theta = thetas(thetaInd);
		fineGrid = linspace(edges(1), edges(end), length(edges)*2);
		[X, Y] = meshgrid(fineGrid);
		R = [cos(theta) -sin(theta); 
		  sin(theta) cos(theta)]
		v = [X(:)'; Y(:)'];
		vRot = R*v;
		Xrot = reshape(vRot(1,:), size(X));
	  end
	  g =1;
	end
  end
  
end