classdef MapComposition < Map
	properties
		Hs
		
	end
	methods
		function this = MapComposition(Hs)
			this.Hs = Hs;
		end
		
		function y = apply_(this, x)
			y = Hs
		end
		
		
		
		function G = makeComposition_(this, H)
			G = MapComposition( [this.Hs, {H}] );
		end
	end
	
	
end