classdef (Abstract) GeneratingFunc
% 2D/3D generating function for x-ray reconstruction

properties
	spaceRadius;
	freqRadius;
	D; % number of dimensions

	separableGroups = []; % e.g. fully separable: [1 2 3], partially separable: [1 1 2], not separable: [] or [1 1 1]
	isotropicGroups = []; % e.g. fully isotropic: [1 1 1], partially isotropic: [1 1 2], not isotrpic: [] or [1 2 3]

end

methods(Abstract)
	phi = val(x)
	Pphi = xray(y, P)
end

methods
	function phiHat = hat(w)
		error('no Fourier transform implemented')
	end
	

	
	
	function checkInput(obj, x)
		if size(x,1) == obj.D
			return;
		elseif size(x,1) == length(unique(obj.isotropicGroups))
			return;
		else
			error('wrong input size, should be %d by N or %d by N', obj.D, length(obj.isotropicDims));
		end
	end
	
	function f = reconstruct(obj, xStep, c)
		phiStart = repmat( -obj.spaceRadius, 1, obj.D);
		phiSize = 2*ceil(obj.spaceRadius ./ xStep )+1;
		
		phiGrid = makeGrid(phiStart, xStep, phiSize); 
		phiVals = obj.val(phiGrid);
		h = reshape(phiVals, phiSize);
		f = imfilter(c, h);
		
	end
	
	
	
	function ySeparableGroups = getProjSeparability(obj, P)
	  % determine the separability of the x-ray transform of phi defined by
	  % the (D-1)xD projection matrix P
	  numGroups = length(unique(obj.separableGroups));
		
		
	  isFunctionOf = P ~= 0; % xi is a function of which yj's?
	  
		% horribly overcomplicated code to determine the merging of the
		% original separable groups...
		groupMembership = false(obj.D-1, numGroups);
		A = zeros(numGroups);
		for d = 1:obj.D-1
			groupMembership(d, obj.separableGroups(isFunctionOf(d,:))) = true;
			for groupInd = find(groupMembership(d,:));
				A(groupInd, groupMembership(d,:)) = 1;
			end
		end
		[V , L] = eigs(diag(sum(A)) - A);
		groups = V(:, sum(L)<eps) ~= 0; % eigenvectors of graph laplacian matching values of 1 give connected comps
		% columns index the new groups after closure, rows are the old groups.

		ySeparableGroups = zeros(1, obj.D-1);
		for d = 1:obj.D-1
			ySeparableGroups(d) = find( groups( find(groupMembership(d,:), 1), : ) );
			
		end
		
		% assign new group numbers starting at 1
		[~, ~, ySeparableGroups(:)] = unique(ySeparableGroups, 'stable');
		
	end
	
	function uniqueProjInds = getUniqueProjs(obj, Ps)
		% determine which x-ray projections defined by the 
		% (D-1) x D projection matrices in Ps are identical due to the isotropy
		% of the kernel
		
		% instead of horribly complicated code, we'll assume 2 or 3D geometry
		
		if length(unique(obj.isotropicGroups)) == 1 % fully isotropic kernel...
			uniqueProjInds = ones(1, size(Ps,3)); % means fully isotropic x-ray proj
			
		elseif length(unique(obj.isotropicGroups)) == 2 && obj.D == 3 % isotropic in one plane 
			groups = unique(obj.isotropicGroups);
			if sum(obj.isotropicGroups == groups(1)) == 2
				nonIsoGroup = groups(2);
			else
				nonIsoGroup = groups(1);
			end
			nonIsoVects = squeeze(Ps(:, obj.isotropicGroups == nonIsoGroup, :));
			[uniqueVects, ~, uniqueProjInds] = unique(nonIsoVects', 'rows');
			
		else % not isotropic at all
			uniqueProjInds = 1:size(Ps,3);
		end
		
	end
	
	
	function S = saveobj(obj)
		S = struct();
		for field = fields(obj)'
			f = field{1};
			S.(f) = obj.(f);
		end
		S.class = class(obj);
	end
	
	
	
end

methods (Static) %% tests -------------------------------------------------
	function test_getProjSeparability()
		phi = BlankTest();
		phi.D = 5;
		phi.separableGroups = [1 2 1 3 3];
		P =[...
			1 0 1 0 0;
			0 0 1 0 0;
			0 0 0 0 1;
			0 1 0 1 0];
		assert(all(phi.getProjSeparability(P) == [1 1 2 2]));
		
		P =[...
			1 1 0 0 0;
			0 0 1 1 0;
			0 0 0 0 1;
			0 1 0 1 0];
		assert(all(phi.getProjSeparability(P) == [1 1 1 1]));
		
		
		phi.separableGroups = [4 5 1 3 2];
		P =[...
			1 0 1 0 0;
			0 0 1 1 0;
			0 0 0 1 0;
			0 1 0 0 0];
		assert(all(phi.getProjSeparability(P) == [1 1 1 2]));
	end
	
	function test_getProjIsotropy()
		phi = BlankTest();
		phi.D = 3;
		phi.isotropicGroups = [1 1 2];
		
		% Ps(:,:,i) is the projection matrix from x-space to y-space (2x3)
		thetas = linspace(0, pi, 10);
		Ps = zeros(2, 3, length(thetas));
		Ps(2, 3, :) = 1;
		for tInd = 1:length(thetas)
			Ps(1, 1, tInd) = cos(thetas(tInd));
			Ps(1, 2, tInd) = sin(thetas(tInd));
		end
		
		phi.getProjIsotropy(Ps);
	end
	
	function obj = loadobj(S)
		const = str2func(S.class);
		obj = const();

		S = rmfield(S, 'class');
		for field = fields(S)'
			f = field{1};
			obj.(f) = S.(f);
		end
		
		
	end
	
end
end

