classdef TestCvgADMM < TestCvg
    % Test convergence by monitoring the primal and dual rediduals in ADMM
    % as described in [1].
    %
    % :param eps_abs: absolute tolerance (default 0, see [1])
    % :param eps_rel: relative tolerance (default 1e-3 see [1])
    %
    % **Reference**
    %
    % [1] Boyd, Stephen, et al. "Distributed optimization and statistical learning via the alternating direction
    % method of multipliers." Foundations and Trends in Machine Learning, 2011.
    %
    % **Warning** the termination criterion described in [1] requires to
    % apply the adjoint of Hn at every iteration, which may be costly
    %
    % **Example** CvOpti=TestCvgADMM(eps_abs,eps_rel)
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
       eps_abs;     % Termination criterion tolerances.
       eps_rel;  
    end
    properties
        evolResPrim
        evolResDual
        count=0;
        yold=[];
    end
    
    methods
        %% Constructor
        function this=TestCvgADMM(eps_abs,eps_rel)
            this.name = 'TestCvgADMM';
            if nargin <2, eps_rel=1e-3; end
            if nargin <1, eps_abs=0; end
            assert(isscalar(eps_abs),'eps_abs must be scalar');
            assert(isscalar(eps_rel),'eps_rel must be scalar');
            this.eps_abs = eps_abs;
            this.eps_rel = eps_rel;
        end
        %% Update method
        function stop = testConvergence(this,opti)
            % Reimplemented from parent class :class:`TestCvg`.
            this.count=this.count+1;
            stop = false;
            
            p = 0; % Number of constraints
            r_norm = 0; % Primal residual norm
            s_norm = 0; % Dual residual norm
            Hnx_norm = 0; % ||Hx||
            y_norm = 0; % ||y||
            adjHnwn_norm = 0; % ||H'w||
            for n = 1:length(opti.wn)
                p = p + length(opti.yn{n});
                r_norm = r_norm + norm(opti.Hnx{n}(:)-opti.yn{n}(:))^2;
                if this.count>1
                    tmp=opti.rho_n(n)*opti.Hn{n}.applyAdjoint(opti.yn{n} - this.yold{n});
                    s_norm = s_norm + norm(tmp(:))^2;
                end
                Hnx_norm = Hnx_norm + norm(opti.Hnx{n}(:))^2;                    
                y_norm = y_norm + norm(opti.yn{n}(:))^2;
                tmp=opti.Hn{n}.applyAdjoint(opti.wn{n});
                adjHnwn_norm = adjHnwn_norm + norm(tmp(:))^2;
            end
            eps_primal = sqrt(p)*this.eps_abs + this.eps_rel*sqrt(max(Hnx_norm, y_norm));
            eps_dual = sqrt(length(opti.xopt))*this.eps_abs + this.eps_rel*sqrt(adjHnwn_norm);
            
            this.evolResPrim(this.count)=gather(r_norm);
            this.evolResDual(this.count)=gather(s_norm);
            
            
            if (this.count>1)&&(sqrt(r_norm) <= eps_primal) && (sqrt(s_norm) <= eps_dual)
                stop = true;
                endingMessage = [this.name,': ADMM convergence criteria reached:',newline, ... 
                    'primal residuals: ' ,num2str(sqrt(r_norm)),' <= ',num2str(eps_primal), newline,... 
                    'dual residuals: ' ,num2str(sqrt(s_norm)),' <= ',num2str(eps_dual)];
                opti.endingMessage = endingMessage;
            else
                this.yold = opti.yn;
            end
            
            
            
        end
    end
end