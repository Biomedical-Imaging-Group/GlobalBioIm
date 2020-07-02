classdef OptiConjGrad < Opti
    % Conjugate gradient optimization algorithm which solves the linear
    % system \\(\\mathrm{Ax=b}\\) by minimizing
    % $$ C(\\mathrm{x})= \\frac12 \\mathrm{x^TAx - b^Tx} $$
    %
    % :param A: symmetric definite positive :class:`LinOp`
    % :param b: right-hand term
    %
    % All attributes of parent class :class:`Opti` are inherited.
    %
    % **Example** CG=OptiConjGrad(A,b,OutOp)
    %
    % See also :class:`Opti`, :class:`OutputOpti` :class:`Cost`

    %%    Copyright (C) 2015
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
    
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        A;  % Linear operator
        b;  % right hand side term
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        r; % residual
        rho_prec;
        p;
    end
    
    methods
        %% Constructor
        function this=OptiConjGrad(A,b)
            this.name='Opti Conjugate Gradient';
            this.A=A;
            assert(checkSize(b,this.A.sizeout),'A sizeout and size of b must be equal');
            this.b=b;
            this.cost=CostL2(this.A.sizeout,0., this.A) - CostLinear(this.A.sizeout, this.b);
        end
        %% Set data b
        function set_b(this,b)
            % Set the right-hand side \\(\\mathrm{b}\\)
            assert(checkSize(b,this.A.sizeout),'A sizeout and size of b must be equal');
            this.b=b;
            
            this.cost.mapsCell{1}.set_y(b);
            this.cost.mapsCell{2}.set_y(b);
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                this.r = this.b - this.A.apply(this.xopt);
                this.p = this.r;
            end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For a detailled
            % algorithm scheme see `here
            % <https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm>`_
            
            rho = dot(this.r(:),this.r(:));
            if this.niter > 0
                beta = rho/this.rho_prec;
                this.p = this.r + beta*this.p;
            end
            q = this.A*this.p;
            alpha = rho/dot(this.p(:), q(:));
            
            % stop if rounding errors
            if dot(this.p(:),this.r(:))<=0
                this.endingMessage = [this.name ,'Rounding errors prevent further optimization'];
                flag = this.OPTI_STOP;
                return
            end
            this.xopt = this.xopt + alpha*this.p;
            this.r = this.r - alpha*q;
            this.rho_prec = rho;
            flag=this.OPTI_NEXT_IT;
        end
    end
end
