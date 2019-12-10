classdef OptiFGP < Opti
    % Fast Gradient Proximal which computes the TV proximity operator which
    % minimizes :class:`Cost` of the form
    % $$ C(\\mathrm{x}) = \\frac12\\|\\mathrm{x} - \\mathrm{y}\\|^2_2 + \\lambda \\|\\mathrm{x} \\|_{TV} $$
    %
    % :param F_0: :class:`CostL2` object
    % :param TV: :class:`CostTV`
    % :param bounds: bounds for set constraint
    % :param gam: descent step (default 1/8)
    % :param lambda: regularization parameter for TV
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
    % **References**
    %
    % [1] Beck, A., and Teboulle, M. (2009). Fast gradient-based algorithms for constrained total variation image denoising 
    % and deblurring problems. IEEE Transactions on Image Processing, 18(11), 2419-2434.
    %
    % **Example** FGP=OptiFGP(F0,TV,bounds)
    %
    % See also :class:`Opti`, :class:`OutputOpti` :class:`Cost`
       
    %% GUI-Header
    % GUInotation-
    % GUIcall-OptiFGP(F0,TV,bounds,gam,lambda)-
    % GUIparam-bounds-vecInt-[]-bounds for set constraint
    % GUIparam-gam-vecInt-[]-descent step (default 1/8)
    % GUIparam-lambda-vecInt-[]-regularization parameter for TV

	%%    Copyright (C) 2015 
    %     T. Pham thanh-an.pham@epfl.ch
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
    properties
        F0; % L2 norm
        TV; % CostTV
        C; % Indicator function for set constraint
        gam; % descent step
    end
    properties (Access = private)
        t;
        F;
        lambda;
        D; % Gradient operator
        ndims;       
        P;
    end
    methods
        function this = OptiFGP(F0,TV,bounds)
            this.name = 'OptiFGP';
            
            assert(isa(F0,'CostL2'), 'F0 must be a CostL2');
            assert(isa(TV,'CostTV') || (isa(TV,'CostMultiplication') && isa(TV.cost2,'CostTV')),...
                'TV must be a CostTV or a CostMultiplication with cost2 as a CostTV');
            if nargin<=2 || isempty(bounds), bounds=[-inf,inf];end
            
            this.F0 = F0;
            if isa(TV,'CostTV')
                this.D = TV.H2;%circular boundary
            else
                this.D = TV.cost2.H2;
            end
            this.ndims = this.D.sizeout(end);
            ElemRep = repmat({':'}, 1, ndims(bounds) - 1 - isvector(bounds));
            this.C = CostRectangle(this.F0.sizein,bounds(ElemRep{:},1),bounds(ElemRep{:},2));
            
            this.TV = TV;
            this.cost = this.F0 + this.TV;
            this.gam = 1/this.D.norm^2;
        end
        
        function setLambda(this,new_l)
            % Set the regularization parameter lambda
            
            if isa(this.TV,'CostTV')
                this.TV = new_l*this.TV;
            else
                this.TV = new_l*this.TV.cost2;
            end
            this.cost = this.F0 + this.TV;
        end
        
        function setBounds(this,new_b)
            % Set constraints bounds
            
            this.bounds = new_b;
            this.C = CostRectangle(this.F0.sizein,this.bounds(1),this.bounds(2));
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                this.xopt = x0;
            end
            if isempty(this.P)
                this.P = zeros_(this.D.sizeout);%dual variable
            end
            this.t = 1;
            
            this.F = this.P;
            if isa(this.TV,'CostTV')
                this.lambda = 1;
            else
                this.lambda = this.TV.cost1;
            end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1].
            
            Pnew = this.F + (this.gam/(this.lambda))*this.D*(this.C.applyProx(this.F0.y - this.lambda*this.D'*(this.F),0));
            Pnew = Pnew./repmat(max(1,sqrt(sum(Pnew.^2, this.ndims + 1))),[ones(1,this.ndims), this.ndims]);%Project L2 ball
            
            tnew = (1 + sqrt(1 + 4*this.t^2))/2;
            this.F = Pnew + (this.t - 1)/tnew*(Pnew - this.P);
            this.t = tnew;
            this.P = Pnew;
            this.xopt = this.C.applyProx(this.F0.y - this.lambda*this.D'*(this.P),0);
            
            flag = this.OPTI_NEXT_IT;
        end
    end
end

