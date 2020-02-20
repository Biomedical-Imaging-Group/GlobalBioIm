classdef MapInversion < Map
    % MapInversion : Builds the inverse :class:`Map`
    %
    % :param M: :class:`Map` object
    %
    % **Example** Minv=MapInversion(M)
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2017 
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
    
    properties (SetAccess = protected,GetAccess = public)
        M;     % Map
    end
    
    %% Constructor
    methods 
        function this = MapInversion(M)
            this.name ='MapInversion';                    
            assert(isa(M,'Map'),'Input should be a Map');
            assert(M.isInvertible,'Input should be invertible');
            this.M = M;
            this.isDifferentiable= this.M.isDifferentiable;
            this.isInvertible=this.M.isInvertible;
            this.sizein =this.M.sizeout;
            this.sizeout =this.M.sizein;   
            this.norm = 1/this.M.norm; 
          end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyInverse_(this,x)
    % - mpower_(this,p)
    % - makeComposition_(this, G)
    methods (Access = protected)
        function y = apply_(this,x) % 
            % Reimplemented from :class:`Map`
            y =this.M.applyInverse(x);
        end
        function y = applyInverse_(this,x)
            % Reimplemented from :class:`Map`
            y =this.M.apply(x);
        end
        function M = mpower_(this,p)
           % Reimplemented from :class:`Map`
            if p==-1
                M=this.M;
            else
                M=mpower_@Map(this,p);
            end
        end
        function M = makeComposition_(this, G)
            % Reimplemented from parent class :class:`Map`.
            if  isequal(this.M,G)
                M=LinOpIdentity(this.sizein);
            else
                M=makeComposition_@Map(this,G);
            end
        end
    end
      
    methods (Access = protected)
        %% Copy
      function this = copyElement(obj)
          this = copyElement@Map(obj);
          this.M = copy(obj.M);
      end
    end
end

