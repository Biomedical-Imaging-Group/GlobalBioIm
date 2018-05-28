classdef OutputOptiConjGrad < OutputOpti
    % OutputOptiADMM class for OptiConjGrad displayings and savings
    %
    % Special :class:`OutpuOpti` for the conjugate gradient algorithm which evaluates
    % $$ C(\\mathrm{x})= \\frac12 \\mathrm{x^TAx - b^Tx} $$
    %
    % See also :class:`Opti` :class:`OutputOpti`

    %%    Copyright (C) 2018
    %     F. Soulez ferreol.soulez@epfl.ch
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
            yty;
        end
            
    methods
        function this=OutputOptiConjGrad(computecost,yty,xtrue,iterVerb)
            if nargin>=1
                if isscalar(computecost)
                    computecost = (computecost ~= 0);
                end
                
                assert(islogical(computecost),'Parameter computecost must be logical');
                this.computecost=computecost;
                this.yty = 0.5.*yty;
            end
            if nargin>=3, this.xtrue=xtrue;end
            if nargin>=4
                assert(isscalar(iterVerb) && iterVerb>=0,'Parameter iterVerb must be a positive integer');
                this.iterVerb=iterVerb;
            end
            if ~isempty(this.xtrue)
                this.isgt=true;
                this.xtrue=xtrue;
                this.normXtrue=norm(this.xtrue(:));
            end
        end
        function cc=computeCost(this,opti)
            % Evaluate the cost function at the current iterate xopt of
            % the given :class:`Opti` opti object  
            r = 0.5.* opti.A.apply(opti.xopt) - opti.b ;
            cc = sum(opti.xopt(:).*r(:)) + this.yty;
        end
    end
end
