classdef OpEWAbs < Map
    % Element-wise Modulus operator: compute the pointwise complex modulus
    % 
    % :param sz: input size
    %
    % All attributes of parent class :class:`Map` are inherited. 
    %
    % **Example** Mod = OpEWAbs(sz)
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2018 
    %     Created: 04/25/2018 (mm/dd/yyyy)
    %     Anthony Berdeu (Laboratoire Hubert Curien)
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
        function this = OpEWAbs(sz)
            this.name ='OpEWAbs ';
            this.sizein=sz;
            this.sizeout=sz;
            this.isDifferentiable=true; % almost ...
            this.isInvertible=false;
        end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(~,x)
            % Reimplemented from parent class :class:`Map`.
            y=abs(x);
        end	
        function x = applyJacobianT_(~,y,v)
            % Reimplemented from parent class :class:`Map`.
            assert(all(v(:)),'Input vector contains zeros');
            x=v./abs(v).*y;
        end	
    end
end