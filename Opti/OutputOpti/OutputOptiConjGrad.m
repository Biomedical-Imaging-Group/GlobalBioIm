classdef OutputOptiConjGrad < OutputOptiSNR
    % OutputOptiConjGrad class displayings and savings dedicated to for OptiConjGrad 
    %
    % The conjugate gradient algorithm minimizes the function
    % $$ C(\\mathrm{x})= \\frac12 \\mathrm{x^TAx - b^Tx} $$
    % However in many cases, it is often used to minimize:
    % $$ F(\\mathrm{x})= \\frac12 \\|H x - y\\|^2_W $$
    % by setting:
    % $$\\mathrm{A} = \\mathrm{H^T W H} \\quad \\text{and}\\quad \\mathrm{b = H^T W y}$$
    % An OutputOptiConjGrad object compute the cost F instead of the cost C. 
    %
    % :param computecost:  boolean, if true the cost function will be computed
    % :param xtrue: ground truth to compute the error with the solution (if provided)
    % :param iterVerb:  message will be displayed every iterVerb iterations (must be a multiple of the :attr:`ItUpOut` parameter of classes :class:`Opti`)
    % :param ytWy:  weighted norm of $y$ : $ \\mathrm{ytWy} = \\mathrm{ y^T\\,W\\,y}$
    %
    % See also :class:`OptiConjGrad` :class:`OutputOpti`
    %

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
            ytWy;
        end
            
    methods
        function this=OutputOptiConjGrad(computecost,yty,xtrue,iterVerb)
            this@OutputOptiSNR(computecost,xtrue,iterVerb);
            if nargin>=1
                this.ytWy = 0.5.*yty;
            end                
        end
        function cc=computeCost(this,opti)
            % Evaluate the cost function at the current iterate xopt of
            % the given :class:`Opti` opti object  
            r = 0.5.* opti.A.apply(opti.xopt) - opti.b ;
            cc = sum(opti.xopt(:).*r(:)) + this.ytWy;
        end
    end
end
