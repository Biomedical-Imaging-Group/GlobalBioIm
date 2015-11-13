classdef KaiserBesselWindow_DPCI < KaiserBesselWindow
	properties
		derivativeDim;
	end
	
	methods
		function obj = KaiserBesselWindow_DPCI(varargin)
			obj = obj@KaiserBesselWindow(varargin{1:end-1});
			if ~isempty(varargin)
				[m, alpha, a, D, dimDim] = deal(varargin{:});
				obj.derivativeDim = dimDim;
			end
		end
		
			
		function p = xray(obj, y, P)
			% projection of generalized kaiser bessel window is now the
			% derivative of the original projection, since this is DPCI
			%
			% "Multidimensional digital image representations using generalized
			% Kaiser-Bessel window functions" (A11)
			
			m1 = obj.m + .5;
			
			ySumSq = sum( y.^2, 1);
			
			isNonzero = ySumSq <= obj.a ^2;
		
			z = obj.alpha * sqrt( 1 - ySumSq / obj.a^2 );
			
			A = -(1 / obj.a) / (obj.alpha^(m1-2) * besseli(m1, obj.alpha));
			
			c1 = zeros(size(ySumSq));
			c1( isNonzero) = y(obj.derivativeDim, isNonzero) / obj.a .* z(isNonzero).^(m1-1) .*besseli(m1-1, z(isNonzero));
			
			% scaling due to projection
			B = obj.a / besseli(obj.m, obj.alpha) * sqrt(2*pi / obj.alpha) * besseli(obj.m + .5, obj.alpha) ;
			
			p = A * B * c1;
			

		end
		
		
	end
	
end