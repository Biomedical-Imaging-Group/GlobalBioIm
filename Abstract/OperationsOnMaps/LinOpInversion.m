classdef LinOpInversion < MapInversion & LinOp
    % LinOpInversion : Builds the inverse :class:`LinOp`
    %
    % :param M: :class:`LinOp` object
    %
    % **Example** Minv=LinOpInversion(M)
    %
    % See also :class:`Map`, :class:`LinOp`, :class:`MapInversion`
    
    %%    Copyright (C) 2017
    %     F. Soulez ferreol.soulez@epfl.ch
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

    %% Constructor
    methods 
        function this = LinOpInversion(M)
            this@MapInversion(M); 
            this.name ='LinOpInversion'; 
            assert(isa(M,'LinOp'),'Input should be a  LinOp');            
          end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyAdjoint_(this,x) 
    % - applyAdjointInverse_
    methods (Access = protected)
        function y = applyAdjoint_(this,x)   
            % Reimplemented from :class:`LinOp`
            y =this.M.applyAdjointInverse(x);
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from :class:`LinOp`
            y =this.M.applyAdjoint(x); 
        end
        % the two function reimplementations below are needed because of
        % the multiple inheritance to specifies which methods to use from
        % parent classes
        function M = mpower_(this,p)
            % Reimplemented from :class:`MapInversion`
            M = mpower_@MapInversion(this,p);
        end
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`MapInversion`
            M = makeComposition_@MapInversion(this,G);
        end
        %% Copy
      	function this = copyElement(obj)
            this = copyElement@LinOp(obj);
            this.M = copy(obj.M);
      	end
    end
end

