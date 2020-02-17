classdef MapComposition < Map
    % MapComposition : Composition of Maps
    % $$ \\mathrm{H}(\\mathrm{x}) = \\mathrm{H}_1 \\left( \\mathrm{H}_2(\\mathrm{x}) \\right) $$
    %
    % :param H1:  left hand side :class:`Map` 
    % :param H2:  right hand side :class:`Map`
    %
    % **Example** H=MapComposition(H1,H2)
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2017
    %     M. McCann michael.mccann@epfl.ch
    %     E. Soubies emmanuel.soubies@epfl.ch
    %
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
    
    properties(SetAccess = protected,GetAccess = public)
        H1;           % Left hand side Map
        H2;           % Right hand side Map
    end
    %% Constructor
    methods
        function this = MapComposition(H1,H2)
            assert(isa(H1,'Map') && isa(H2,'Map'),'H1 and H2 have to be a Map');
            this.H1 = H1;
            this.H2 = H2;
            % Set Properties
            % norm
            if this.H1.norm ~= -1 && this.H2.norm ~= -1
                this.norm = this.H1.norm * this.H2.norm;
            else
                this.norm=-1;
            end
            % isInvertible
            if this.H1.isInvertible && this.H2.isInvertible, this.isInvertible=true;end
            % isDifferentiable
            if this.H1.isDifferentiable && this.H2.isDifferentiable, this.isDifferentiable=true; end
            % sizein/out
            this.sizein=H2.sizein;
            this.sizeout=H1.sizeout;
            this.name=sprintf('MapComposition( %s ; %s )',H1.name,H2.name);      
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - applyInverse_(this,y)
	methods (Access = protected)
        function y = apply_(this, x)
            % Reimplemented from :class:`Map`
            y = this.H1.apply(this.H2.apply(x));
        end
        function x = applyJacobianT_(this, y, v)
            % Reimplemented from :class:`Map`
            if this.isDifferentiable
                x=this.H2.applyJacobianT(this.H1.applyJacobianT(y,this.H2.apply(v)),v);
            else
                x = applyJacobianT_@Map(this,y,v);
            end
        end
        function x = applyInverse_(this, y)
            % Reimplemented from :class:`Map`
            if this.isInvertible
                x=this.H2.applyInverse(this.H1.applyInverse(y));
            else
                x = applyInverse_@Map(this,y);
            end
        end
        function M = makeComposition_(this,G)
             % Reimplemented from :class:`Map`
			 H2G = this.H2*G; % Since H1*H2 is not simplified, first try to compose H2 with G
			 if ~isa(H2G, 'MapComposition')
				 M=this.H1*H2G;
			 else
				 M = MapComposition(this, G);
			 end
        end
    end	
    
    methods (Access = protected)
        %% Copy
      function this = copyElement(obj)
          this = copyElement@Map(obj);
          this.H1 = copyElement@Map(obj.H1);
          this.H2 = copyElement@Map(obj.H2);
      end
    end
end