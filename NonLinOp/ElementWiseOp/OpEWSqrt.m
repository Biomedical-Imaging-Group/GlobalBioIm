classdef OpEWSqrt < Map
    % Element-wise square root operator: 
    % $$[\\mathrm{H}(x)]_k = \\sqrt{(x_k)} $$
    % 
    % :param sz: input size
    %
    % All attributes of parent class :class:`Map` are inherited. 
    %
    % **Example** H=OpEWSquareRoot(sz)
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
    
    %% Constructor
    methods
        function this = OpEWSqrt(sz)
            this.name ='OpEWSqrt';
            this.sizein=sz;
            this.sizeout=sz;
            this.isDifferentiable=true;
            this.isInvertible=true;
        end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`Map`.

            assert(sum(x(:)<0)==0,'Input vector contains negative values');
            y=sqrt(x);
        end	
        function x = applyJacobianT_(this,y,v)
            % Reimplemented from parent class :class:`Map`.

            x=1./(2*sqrt(v)).*y;
        end	
    end
end