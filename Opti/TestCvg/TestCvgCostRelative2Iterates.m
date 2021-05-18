classdef TestCvgCostRelative2Iterates  < TestCvg
    % TestCvgCostRelative2Iterates stops the optimization when the relative decrease of cost with respect to the last 2 iterations is below the value COSTRELATIVETOL
    %
    % :param costRelativeTol:  relative tolerance on cost function
    % :param costIndex: select a specific cost function among a sum in the case where the optimized cost function is a sum of cost functions
    %
    % **Example** CvOpti=TestCvgCostRelative2Iterates(costRelativeTol, costIndex )
    %
    % See also :class:`TestCvg`
    
    %%    Copyright (C) 2021
    %     T. Debarre thomas.debarre@epfl.ch
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
        costRelativeTol=1e-5;      % stopping criteria tolerance on the relative difference btw the current iterate and the 2 previous ones
    end
    properties (SetAccess = protected,GetAccess = protected)
        oldCost=[];               % cost from previous iterate
        oldOldCost=[];            % cost from 2 iterates ago 
        costIndex=0;              % index of the cost function
    end
    methods
        %% Constructor
        function this=TestCvgCostRelative2Iterates( costRelativeTol, costIndex)
            this.name = 'TestCvgCostRelative2Iterates';
            assert(isscalar(costRelativeTol),'costRelativeTol must be scalar');
            this.costRelativeTol =costRelativeTol;
            if(nargin==2)
                %                assert(isscalar(costIndex) ,'costIndex must be a scalar integer');
                this.costIndex = costIndex;
            end
        end
        %% Update method
        function stop = testConvergence(this,opti)
            % Tests algorithm convergence from the relative difference
            % between the current cost and that of the two previous
            % iterations. Note: this criterion can be useful for algorithms
            % which can have 2 equal consecutive iterates before
            % convergence, such as FISTA (whose updates depend on the
            % previous 2 iterates).
            %
            % :return: boolean true if
            %     $$ \max \left( \\frac{\\left| C(\\mathrm{x}^{k}) - C(\\mathrm{x}^{k-1})\\right|}{\\left|C(\\mathrm{x}^{k-1})\\right|}, \\frac{\\left| C(\\mathrm{x}^{k}) - C(\\mathrm{x}^{k-2})\\right|}{\\left|C(\\mathrm{x}^{k-2})\\right|} \right) < \\mathrm{costRelativeTol}$$
            
            
            stop = false;
            if ( isa(opti.cost,'CostSummation')&& all(this.costIndex>0) &&all(this.costIndex<=opti.cost.numMaps) )
                f = 0;
                for n=1:numel(this.costIndex)
                    f = f+opti.cost.mapsCell{this.costIndex(n)}*opti.xopt;
                end
            else
                f = opti.cost*opti.xopt;
            end
            
            if ~isempty(this.oldCost) && ~isempty(this.oldOldCost)
                r = max(abs( this.oldCost - f)./this.oldCost, abs(this.oldOldCost - f)./this.oldOldCost);
                if( r < this.costRelativeTol)
                    stop = true;
                    endingMessage = [this.name,': Cost variation between current iterate and the two previous ones below the relative tolerance : ',num2str(r),' < ',num2str(this.costRelativeTol)];
                    opti.endingMessage = endingMessage;
                end
            end
            this.oldOldCost = this.oldCost;
            this.oldCost = f;
            
        end
    end
end
