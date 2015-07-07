classdef Prox < handle
    %% Prox : Moreau Proximity Operator generic class
    %  Matlab Inverse Problems Library
    % The Prox meta class implement generic methods for all Moreau Proximity 
    % operators:
    %
    % $$ \mathrm{Prox}_{\alpha\,\phi}(\mathbf{u}) = \mathrm{arg\,min}_{\mathbf{x}}\  \alpha\,\phi(\mathbf{x}) +\frac{1}{2} || \mathbf{u} - \mathbf{x}||_{2}^{2} $$
    %
    %% Properties
    % * |name|          - name of the function  $\phi$ 
    % * |alpha|         - hyperparameter $\alpha$
    % * |sz|            - size of x
    %
    %% Methods
    % * |Apply|    - Apply the proximity operator of $\phi$
    % * |Cost|     - Compute the total cost $\phi(\mathrm{Prox}_{\alpha\,\phi}(\mathbf{x})) + \frac{1}{2} || \mathrm{Prox}_{\alpha\,\phi}(\mathbf{x})) - \mathbf{x}||_{2}^{2} $
    % * |FCost|    - Compute $\phi(\mathrm{Prox}_{\alpha\,\phi}(\mathbf{x}))$
    % * |Residual| - Compute the residual  % $\mathrm{Prox}_{\alpha\,\phi}(\mathbf{x}) - \mathbf{x}$
    %%
    
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
        alpha ;         % hyperparameter $\alpha$
        sz;             % dimensions of the vector space
    end
    
    methods
        function Apply(~,~,~) % Apply the prox
            error('Operator not implemented');
        end
        function Cost(~,~) % Get the global cost
            error('Cost not implemented');
        end
        function FCost(~,~) % get the function cost
            error('fCost not implemented');
        end
        function Residual(~,~) % get the residuals
            error('residuals not implemented');
        end
    end
end
