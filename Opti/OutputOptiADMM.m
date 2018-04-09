classdef OutputOptiADMM < OutputOpti
    % OutputOptiADMM class for OptiADMM displayings and savings
    %
    % Special :class:`OutpuOpti` for ADMM which evaluates the Lagrangian
    % instead of the cost function.
    %
    % See also :class:`Opti` :class:`OutputOpti`

    %%    Copyright (C) 2018
    %     E. Soubies emmanuel.soubies@epfl.ch
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
    
    methods
        %% Constructor
        function this=OutputOptiADMM(computecost,xtrue,iterVerb)
            this@OutputOpti(computecost,xtrue,iterVerb);
        end
        %% computeCost  method
        function cc=computeCost(this,opti)
            % Reimplemented from parent class :class:`OutputOpti`.
            % Evaluates the Lagrangian instead of the objective function.
            if isempty(opti.F0)
                cc=0;
            else
                cc = opti.F0*opti.xopt;
            end
            for n=1:numel(opti.Fn)
                tmp=opti.Hn{n}*opti.xopt-opti.yn{n}+opti.wn{n}/opti.rho_n(n);
                cc = cc+opti.Fn{n}*opti.yn{n} + 0.5*opti.rho_n(n)*norm(tmp(:),'fro')^2;
            end
        end
    end
end
