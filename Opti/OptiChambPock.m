classdef OptiChambPock < Opti
    % Chambolle-Pock optimization algorithm [1] which minimizes :class:`Cost` of the form
	% $$ C(\\mathrm{x}) = F(\\mathrm{Hx}) + G(\\mathrm{x}) $$
    %
    % :param F: a :class:`Cost` with an implementation of the :meth:`prox`.
    % :param G: a :class:`Cost` with an implementation of the :meth:`prox`.
    % :param H: a :class:`LinOp`.
    % :param tau: parameter of the algorithm (default 1)
    % :param sig: parameter of the algorithm which is computed automatically if H.norm is different from -1.
    % :param var: select the "bar" variable of the algorithm (see [1]):
    %
    %   - if 1 (default) then the primal variable \\(\\bar{\\mathrm{x}} = 2\\mathrm{x}_n  - \\mathrm{x}_{n-1}\\) is used 
    %
	%   - if 2 then the dual variable \\(\\bar{\\mathrm{y}} = 2\\mathrm{y}_n  - \\mathrm{y}_{n-1} \\) is used
    %
    % :param gam: if non-empty, accelerated version (see [1]). Here, \\(G\\) or \\(F^*\\) is uniformly convex with \\(\\nabla G^*\\) or \\(\\nabla F\\) 1/gam-Lipschitz:
    %
    %   - If \\(G\\) is uniformly convex then set the parameter var to 1. 
    %
    %   - If \\(F^*\\) is uniformly convex then set the parameter var to 2
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
	% **Note-1**: In fact, \\(F\\) needs only the prox of its fenchel transform :meth:`prox_fench` (which is implemented 
	% as soon as $F$ has an implementation of the prox, see :class:`Cost`).
    %
    % **Note-2**:
    %
    %   - To ensure convergence (see [1]), parameters sig and tau have to verify 
    %     $$ \\sigma \\times \\tau \\times \\Vert \\mathrm{H} \\Vert^2 < 1 $$ 
    %     where \\(\\Vert \\mathrm{H}\\Vert\\) denotes the norm of the linear operator H.
    %
	%   - When the accelerated version is used (i.e. parameter gam is non-empty), 
    %     sig and tau will be updated at each iteration and the initial 
    %     ones (given by user)  have to verify
    %     $$ \\sigma \\times \\tau \\times \\Vert \\mathrm{H} \\Vert^2 \\leq 1 $$
    %
    % **Reference**: 
    %
    % [1] Chambolle, Antonin, and Thomas Pock. "A first-order primal-dual algorithm for convex problems with 
	% applications to imaging." Journal of Mathematical Imaging and Vision 40.1, pp 120-145 (2011).	
    %
    % **Example** CP=OptiChambPock(F,H,G)
    %
    % See also :class:`Opti` :class:`OutputOpti` :class:`Cost`    
    
    %%    Copyright (C) 2017 
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

    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		F;  % Cost F
		G;  % Cost G
		H;  % LinOp H
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = public)
		theta=1; % parameter of the algorithm (fixed to 1 but can could be set as a user defined 
                 % parameter for future upgrades of the Library)
        % Internal variables
        y;
        xbar;
        Kxbar;
        Kxopt;
        Kxold;
        yold;
        ybar;
        KTyold;
        KTy;
        KTybar;
    end
    % Full public properties
    properties
		sig=[];   % parameter of the algorithm
		tau=1;    % parameter of the algorithm
		gam=[];   % parameter of the algorithm for accelerated version
		var=1;    % select the bar variable in the algorith (primal 1/ dual 2)
    end
    
    methods
    	%% Constructor
    	function this=OptiChambPock(F,H,G)
    		this.name='Opti Chambolle-Pock';
    		this.cost=F*H+G;
    		this.F=F;
    		this.G=G;
            this.H=H;
            this.needxold =true;
            if this.H.norm>=0
                this.sig=1/(this.tau*this.H.norm^2)-eps;
            end
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~isempty(x0) % To restart from current state if wanted
                assert(~isempty(this.sig),'parameter sig is not setted');
                assert(~isempty(this.tau),'parameter tau is not setted');
                this.y=this.H.apply(x0);
                if this.var==1
                    this.xbar=x0;
                    this.Kxbar= this.y;
                else
                    this.ybar= this.y;
                    this.KTybar=this.H.applyAdjoint(this.ybar);
                    this.KTy=this.KTybar;
                end
                this.Kxopt=this.y;
            end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1].

            if this.var==1 % === using xbar
                this.Kxold=this.Kxopt;
                % - Algorithm iteration
                this.y=this.F.applyProxFench(this.y+this.sig*this.Kxbar,this.sig);
                this.xopt=this.G.applyProx(this.xopt-this.tau*this.H.applyAdjoint(this.y),this.tau);
                this.Kxopt=this.H.apply(this.xopt);
                if ~isempty(this.gam) % acceleration => uodate theta, tau and sig according to [1]
                    this.theta=1/sqrt(1+2*this.gam*this.tau);
                    this.tau=this.theta*this.tau;
                    this.sig=this.sig/this.theta;
                end
                this.xbar=this.xopt+this.theta*(this.xopt - this.xold);
                this.Kxbar=this.Kxopt+this.theta*(this.Kxopt - this.Kxold);
            else % === using ybar
                this.yold=this.y;
                this.KTyold=this.KTy;
                % -- Algorithm iteration
                this.xopt=this.G.prox(this.xopt-this.tau*this.KTybar,this.tau);
                this.Kxopt=this.H.apply(this.xopt);
                this.y=this.F.prox_fench(this.y+this.sig*this.Kxopt,this.sig);
                this.KTy=this.H.adjoint(this.y);
                if ~isempty(this.gam) % acceleration => uodate theta, tau and sig according to [1]
                    this.theta=1/sqrt(1+2*this.gam*this.sig);
                    this.tau=this.tau/this.theta;
                    this.sig=this.sig*this.theta;
                end
                this.ybar=(1+this.theta)*this.y - this.theta*this.yold;
                this.KTybar=(1+this.theta)*this.KTy - this.theta*this.KTyold;
            end
            flag=this.OPTI_NEXT_IT;
        end
	end
end
