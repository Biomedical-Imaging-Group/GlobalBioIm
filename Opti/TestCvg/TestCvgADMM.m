classdef TestCvgADMM < TestCvg
    % Test convergence by monitoring the primal and dual rediduals in ADMM
    %
    % **Reference**
    %
    % [1] Boyd, Stephen, et al. "Distributed optimization and statistical learning via the alternating direction
    % method of multipliers." Foundations and Trends in Machine Learning, 2011.
    %
    % Warning: the termination criterion described in [1] requires to
    % apply the adjoint of Hn at every iteration, which may be costly
    %
    % **Example** CvOpti=TestCvgTemplate()
    %
    % **Important** The update method should have an unique imput that is the :class:`Opti` object in order to
    % be generic for all Optimization routines. Hence the update method has access (in reading mode)
    % to all the properties of :class:`Opti` objects.
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
    
    properties (SetAccess = protected,GetAccess = public)
    end
    properties (SetAccess = public,GetAccess = public)
       eps_abs=0;   % Termination criterion tolerances.
        eps_rel=0;  
    end
    
    methods
        %% Constructor
        function this=TestCvgADMM(eps_abs,eps_rel)
            this.name = 'TestCvgADMM';
assert(isscalar(eps_abs),'eps_abs must be scalar');
assert(isscalar(eps_rel),'eps_rel must be scalar');
this.eps_abs = eps_abs;
this.eps_rel = eps_rel;
        end
        %% Update method
        function stop = testConvergence(this,opti)       
            stop = false;
         
                p = 0; % Number of constraints
                r_norm = 0; % Primal residual norm
                s_norm = 0; % Dual residual norm
                Hnx_norm = 0; % ||Hx||
                y_norm = 0; % ||y||
                adjHnwn_norm = 0; % ||H'w||
                for n = 1:length(opti.wn)
                    p = p + length(opti.yn{n});
                    r_norm = r_norm + norm(opti.Hnx{n}-opti.yn{n})^2;
                    s_norm = s_norm + norm(opti.rho_n(n)*opti.Hn{n}.applyAdjoint(opti.yn{n} - opti.yold{n}))^2;
                    Hnx_norm = Hnx_norm + norm(opti.Hnx{n})^2;
                    y_norm = y_norm + norm(opti.yn{n})^2;
                    adjHnwn_norm = adjHnwn_norm + norm(opti.Hn{n}.applyAdjoint(opti.wn{n}))^2;
                end
                eps_primal = sqrt(p)*this.eps_abs + this.eps_rel*sqrt(max(Hnx_norm, y_norm));
                eps_dual = sqrt(length(opti.xopt))*this.eps_abs + this.eps_rel*sqrt(adjHnwn_norm);
                
            if (sqrt(r_norm) <= eps_primal) && (sqrt(s_norm) <= eps_dual)
                stop = true;
                endingMessage = [this.name,' convergence reached'];
                opti.endingMessage = endingMessage;
            end
        end
    end
end