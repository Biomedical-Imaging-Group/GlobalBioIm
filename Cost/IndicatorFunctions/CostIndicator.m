classdef CostIndicator < Cost
    %% CostIndicator : Indicator function
    %  Matlab Inverse Problems Library
    % The CostIndicator meta class implements generic methods for all
    % indicator cost functions
    %
    %
    %     Copyright (C) 2017 F. Soulez ferreol.soulez@epfl.ch
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
    end
    
    
    methods
        %% Gradient of the Functional
        function grad(~,~)
            error('Cannot evaluate gradient of an indicator function');
        end
        %% Compute the functional and its gradient
        function  eval_grad(~,~)       
            error('Cannot evaluate gradient of an indicator function');
        end
    end
end
