classdef TestCvgCostRelative  < TestCvg
    % TestCvgCostRelative stops the optimization when the cost function is below the value COSTRELATIVETOL
    %
    % :param costRelativeTol:  relative tolerance on cost function
    % :param costIndex: select a specific cost function among a sum in the case where the optimized cost function is a sum of cost functions
    %
    % **Example** CvOpti=TestCvgCostRelative(costRelativeTol, costIndex )
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
        costRelativeTol=1e-5;      % stopping criteria tolerance on the relative difference btw two successive iterates
    end
    properties (SetAccess = protected,GetAccess = protected)
        oldCost=[];
        costIndex=0;              % index of the cost function
    end
    methods
        %% Constructor
        function this=TestCvgCostRelative( costRelativeTol, costIndex)
            this.name = 'TestCvgCostRelative';
            assert(isscalar(costRelativeTol),'costRelativeTol must be scalar');
            this.costRelativeTol =costRelativeTol;
            if(nargin==2)
                %                assert(isscalar(costIndex) ,'costIndex must be a scalar integer');
                this.costIndex = costIndex;
            end
        end
        %% Update method
        function stop = testConvergence(this,opti)
            % Tests algorithm convergence from the relative difference between two successive value of the cost function
            %
            % :return: boolean true if
            %     $$ \\frac{\\left| C(\\mathrm{x}^{k}) - C(\\mathrm{x}^{k-1})\\right|}{\\left|C(\\mathrm{x}^{k-1})\\right|} < \\mathrm{costRelativeTol}$$
            
            
            stop = false;
            if ( isa(opti.cost,'CostSummation')&& all(this.costIndex>0) &&all(this.costIndex<=opti.cost.numMaps) )
                f = 0;
                for n=1:numel(this.costIndex)
                    f = f+opti.cost.mapsCell{this.costIndex(n)}*opti.xopt;
                end
            else
                f = opti.cost*opti.xopt;
            end
            
            
            
            if ~isempty(this.oldCost)
                r = abs( this.oldCost - f)./abs(this.oldCost);
                if( r < this.costRelativeTol)
                    stop  =true;
                    endingMessage = [this.name,': Cost variation between two successive iterates below the relative tolerance : ',num2str(r),' < ',num2str(this.costRelativeTol)];
                    opti.endingMessage = endingMessage;
                end
            end
            this.oldCost = f;
        end
    end
end
