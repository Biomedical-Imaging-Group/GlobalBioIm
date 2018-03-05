classdef TestCvgCostAbsolute  < handle
    % TestCvgCostAbsolute stops the optimization when the cost function is below the value COSTABSOLUTETOL
    %
    % :param verbose: if true will display a message before stopping the algorithm.
    % :param costAbsoluteTol:  absolute tolerance on cost function
    % :param costIndex: select a specific cost function among a sum in the case where the optimized cost function is a sum of cost functions
    %
    % **Example** CvOpti=TestCvgCostAbsolute(verbose,costAbsoluteTol, costIndex )
    %
    % See also :class:`Opti`
    
    %%    Copyright (C) 2018
    %     F. Soulez ferreol.soulez@univ-lyon1.fr
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
    
    properties (Constant)
        name = 'TestCvgCostAbsolute'% name of the optimization algorithm
    end
    properties (SetAccess = protected,GetAccess = public)
        verbose = true
        count=1; % counter
    end
    properties (SetAccess = public,GetAccess = public)
        costAbsoluteTol=1e-5;      % stopping criteria tolerance on the value of the cost function
        costIndex=0;              % index of the cost function
    end
    methods
        %% Constructor
        function this=TestCvgCostAbsolute(verbose, costAbsoluteTol, costIndex)
            this.verbose= verbose;
            assert(isscalar(costAbsoluteTol),'costAbsoluteTol must be scalar');
            this.costAbsoluteTol =costAbsoluteTol;
            if(nargin==3)
                assert(isscalar(costIndex) && isinteger(costIndex),'costIndex must be a scalar integer');
                this.costIndex = costIndex;
            end
        end
        %% Update method
        function stop = testConvergence(this,opti)            
            % Tests algorithm convergence 
            %
            % :returns stop: boolean true if Cost <costAbsoluteTol
            stop = false;
            if (this.costIndex && isa(opti.cost,'CostSummation') && this.costIndex<= opti.cost.numMaps)
                f = opti.cost.mapsCell{this.costIndex}*opti.xopt;
            else
                f = opti.cost*opti.xopt;
            end
            if( f < this.costAbsoluteTol)
                stop  =true;
                message = [this.name,': Cost below the absolute tolerance : ',f,' < ',this.costAbsoluteTol];
                this.message = message;
                if this.verbose
                    disp(message)
                end
            end
        end
    end
end
