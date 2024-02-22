classdef OptiDouglasRachford < Opti
	% Douglas Rachford splitting optimization algorithm which minimizes
	% :class:`Cost` of the form
    % $$ C(\\mathrm{x}) = F_1(\\mathrm{x}) + F_2(\\mathrm{L x}) $$
    %
    % :param F_1: a :class:`Cost` with an implementation of the :meth:`prox`.
    % :param F_2: a :class:`Cost` with an implementation of the :meth:`prox`.
    % :param L: a :class:`LinOp` such that \\(\\mathrm{LL^T} = \\nu \\mathrm{I} \\)
	% :param gamma: \\(\\in [0,+\\inf[\\)
	% :param lambda: \\(\\in ]0,2[\\) the relaxation parmeter
    %
    % All attributes of parent class :class:`Opti` are inherited. 
	%
    % **Example** DR=OptiDouglasRachford(F1, F2, L, gamma, lambda)
    %
	% See also :class:`Opti`, :class:`OutputOpti`, :class:`Cost`

    %%    Copyright (C) (2017)
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
		F1;
		F2;
		lambda
		gamma
    end
    properties (SetAccess = protected,GetAccess = protected)
        y;    
    end
    
	methods
		function this = OptiDouglasRachford(F1, F2, L, gamma, lambda)
            if isempty(L), L=LinOpIdentity(F2.sizein); end;
            this.name='OptiDouglasRachford';
			% F1, F2
			this.F1 = F1;
			this.F2 = F2*L;
            this.cost=this.F1+this.F2;
			
			% gamma
			if numel(gamma)==1
				gamma = [gamma, gamma];
			end
			this.gamma = gamma;
			
			% lambda
			this.lambda = lambda;
			
		end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                this.y = x0;
            end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`.

            this.xopt = this.F2.applyProx(this.y, this.gamma(2));
            this.y = this.y + this.lambda .* ( this.F1.applyProx(2.*this.xopt- this.y,this.gamma(1)) - this.xopt);
            flag=this.OPTI_NEXT_IT;
        end
    end
end
