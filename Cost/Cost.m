classdef Cost < handle
    %% Cost : Cost function generic class
    %  Matlab Inverse Problems Library
    % The Cost meta class implement generic methods for all cost functions
    %
    %% Properties
    % * |name|          - name of the function  
    % * |sizeIn|          - size of input space (kernel)
    %
    %% Methods
    % * |update|        - Update the input x
    % * |getCost|       - compute the cost 
    % * |getGradient|   - Compute $\phi(\mathrm{Prox}_{\alpha\,\phi}(\mathbf{x}))$
    % * |prox|          - Compute the  Moreau proximal mapping operator
    %%
    
    
    %     Copyright (C) 2016 F. Soulez ferreol.soulez@epfl.ch
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
        name = 'none'   % name of the linear operator
        sizeIn;             % dimensions of the input vector space
    end
    
    methods
        function update(~,~,~) % 
            error('Operator not implemented');
        end
        function getCost(~,~) % Get the  cost
            error('Cost not implemented');
        end
        function getGradient(~,~) % get the function cost
            error('fCost not implemented');
        end
        function prox(~,~) % Apply the prox
            error('residuals not implemented');
        end
    end
end
