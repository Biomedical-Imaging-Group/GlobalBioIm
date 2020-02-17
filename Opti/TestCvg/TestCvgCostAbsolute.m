classdef TestCvgCostAbsolute  < TestCvg
    % TestCvgCostAbsolute stops the optimization when the cost function is below the value COSTABSOLUTETOL
    %
    % :param costAbsoluteTol:  absolute tolerance on cost function
    % :param costIndex: select a specific cost function among a sum in the case where the optimized cost function is a sum of cost functions
    %
    % **Example** CvOpti=TestCvgCostAbsolute(costAbsoluteTol, costIndex )
    %
    % See also :class:`TestCvg`
    
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
    
    properties (SetAccess = public,GetAccess = public)
        costAbsoluteTol=1e-5;      % stopping criteria tolerance on the value of the cost function
        costIndex=0;              % index of the cost function
    end
    methods
        %% Constructor
        function this=TestCvgCostAbsolute( costAbsoluteTol, costIndex)
            this.name ='TestCvgCostAbsolute';

            assert(isscalar(costAbsoluteTol),'costAbsoluteTol must be scalar');
            this.costAbsoluteTol =costAbsoluteTol;
            if(nargin==2)
                % assert(isscalar(costIndex),'costIndex must be a scalar integer');
                this.costIndex = costIndex;
            end
        end
        %% Update method
        function stop = testConvergence(this,opti)            
            % Tests algorithm convergence 
            %
            % :return: boolean true if \\( C(\\mathrm{x^k}) < \\mathrm{costAbsoluteTol}\\)
            
            stop = false;
            if ( isa(opti.cost,'CostSummation')&& all(this.costIndex>0) &&all(this.costIndex<=opti.cost.numMaps) )
                f = 0;
                for n=1:numel(this.costIndex)
                    f = f+opti.cost.mapsCell{this.costIndex(n)}*opti.xopt;
                end
            else
                f = opti.cost*opti.xopt;
            end

            if( f < this.costAbsoluteTol)
                stop  =true;
                endingMessage = [this.name,': Cost below the absolute tolerance : ',num2str(f),' < ',num2str(this.costAbsoluteTol)];
                opti.endingMessage = endingMessage;
            end
        end
    end
end
