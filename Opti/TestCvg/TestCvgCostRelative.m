classdef TestCvgCostRelative  < handle
    % TestCvgCostRelative stops the optimization when the cost function is below the value COSTABSOLUTETOL
    %
    % :param verbose: if true will display a message before stopping the algorithm.
    % :param costRelativeTol:  absolute tolerance on cost function
    % :param costIndex: select a specific cost function among a sum in the case where the optimized cost function is a sum of cost functions
    %
    % **Example** CvOpti=TestCvgCostRelative(verbose,costRelativeTol, costIndex )
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
    
    properties (SetAccess = protected,GetAccess = public)
        name = 'TestCvgCostRelative'% name of the optimization algorithm
        verbose = true
        count=1; % counter
    end
    properties (SetAccess = public,GetAccess = public)
        costRelativeTol=1e-5;      % stopping criteria tolerance on the relative difference btw two successive iterates
        costIndex=0;              % index of the cost function
    end
    properties (SetAccess = protected,GetAccess = protected)
        oldCost=[];      % stopping criteria tolerance on the relative difference btw two successive iterates
    end
    methods
        %% Constructor
        function this=TestCvgCostRelative(verbose, costRelativeTol, costIndex)
            this.verbose= verbose;
            assert(isscalar(costRelativeTol),'costRelativeTol must be scalar');
            this.costRelativeTol =costRelativeTol;
            if(nargin==3)
                assert(isscalar(costIndex) && isinteger(costIndex),'costIndex must be a scalar integer');
                this.costIndex = costIndex;
            end
        end
        %% Update method
        function stop = testConvergence(this,opti)           
            % Tests algorithm convergence from the relative difference between two successive value of the cost function
            % :returns stop: boolean true if 
            % $$ \\frac{\\| \\mathrm{f}^{k} - \\mathrm{f}^{k-1}\\|}{\\|\\mathrm{f}^{k-1}\\|} < \\text{costRelativeTol}.$$
            
            r=this.xopt-this.xold;
            xdiff=norm()/(norm(this.xold(:))+eps);
            stop=xdiff<this.xtol;
            stop = false;
            if (this.costIndex && isa(opti.cost,'CostSummation') && this.costIndex<= opti.cost.numMaps)
                f = opti.cost.mapsCell{this.costIndex}*opti.xopt;
            else
                f = opti.cost*opti.xopt;
            end
            
            
            
            if ~isempty(this.oldCost)
                r = (f - this.oldCost)./this.oldCost;
                if( r < this.costRelativeTol)
                    stop  =true;
                    message = [this.name,': Cost variation between two successive iterates below the relative tolerance : ',f,' < ',this.costRelativeTol];
                    this.message = message;
                    if this.verbose
                        disp(message)
                    end
                end
            end
            this.oldCost = f;
        end
    end
end
