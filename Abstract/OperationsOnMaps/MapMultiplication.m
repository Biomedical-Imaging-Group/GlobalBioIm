classdef MapMultiplication < Map
    % MapMultiplication: Multiplication of Maps
    % $$ \\mathrm{H}(\\mathrm{x}) =  \\mathrm{H}_1(\\mathrm{x}) \\times \\mathrm{H}_2(\\mathrm{x}) $$
    %
    % :param Map1:  :class:`Map` object
    % :param Map2:  :class:`Map` object
    %
    % **Example** H=MapMultiplication(Map1,Map2)
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2018
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
        M1;
        M2;
    end
    %% Constructor
    methods
        function this = MapMultiplication(Map1,Map2)
            this.name ='MapMultiplication';	
            assert((isa(Map1,'Map') && isa(Map2,'Map')),'Arguments should be Map objects');
            this.M1=Map1;
            this.M2=Map2;
            this.isDifferentiable= Map1.isDifferentiable*Map2.isDifferentiable;
            this.isInvertible=false;
            this.sizein = this.M1.sizein;
            this.sizeout = this.M2.sizeout;
            assert(cmpSize(this.sizein,this.M2.sizein),'Given Maps do not have consistent  sizein') ;
            assert(cmpSize(this.sizeout,this.M2.sizeout),'Given Maps do not have the consistent sizeout ');
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - makeComposition_(this,G)
    methods (Access = protected)
        function y = apply_(this,x) 
            % Reimplemented from :class:`Map`   
            y = (this.M1*x).*(this.M2*x);
        end  
        function x = applyJacobianT_(this, y, v)
            % Reimplemented from :class:`Map`   
            
            x=(this.M1*v).*this.M2.applyJacobianT(y,v) + (this.M2*v).*this.M1.applyJacobianT(y,v);
        end     
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`Map`  
            M=MapMultiplication(this.M1*G,this.M2*G);
        end
    end
    
    
    methods (Access = protected)
        %% Copy
      function this = copyElement(obj)
          this = copyElement@Map(obj);
          this.M1 = copyElement@Map(obj.M1);
          this.M2 = copyElement@Map(obj.M2);
      end
end

