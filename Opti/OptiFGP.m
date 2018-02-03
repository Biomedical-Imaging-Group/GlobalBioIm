classdef OptiFGP < Opti
    
    %OptiFGP optimizes the TV proximity operator using Fast Gradient
    %Proximal.
    properties
        F0;
        TV;
        OpTV;
        C; % Indicator function for set constraint
        D; % Gradient operator
        L; %L: Lipschitz constant of F0 or others
        ndims;
        P;
        Nesterov = 0; %boolean for nesterov acceleration (FGP)
    end
    methods
        function this = OptiFGP(F0,TV,bounds,OutOp)
            this.name = 'OptiFGP';
            
            assert(isa(F0,'CostL2'), 'F0 must be a CostL2');
            assert(isa(TV,'CostTV') || (isa(TV,'CostMultiplication') && isa(TV.cost2,'CostTV')),...
                'TV must be a CostTV or a CostMultiplication with cost2 as a CostTV');
            if nargin<=2 || isempty(bounds), bounds=[-inf,inf];end
            if nargin<=3 || isempty(OutOp), OutOp = OutputOpti(0);end
            
            this.OutOp = OutOp;
            this.F0 = F0;
            if isa(TV,'CostTV')
                this.D = TV.H2;%circular boundary
                TV = 1*TV;
            else
                this.D = TV.cost2.H2;
            end
            this.ndims = this.D.sizeout(end);
            this.C = CostRectangle(this.F0.sizein,bounds(1),bounds(2));
            this.OpTV = LinOpTV(this.F0.sizein,this.D.bc); %For the dual variables
            this.TV = TV;
            this.cost = this.F0 + this.TV;
            if this.F0.lip~=-1
                this.L = 8*this.F0.lip;
            else
                this.L = 8;
            end
        end
        
        function setLambda(this,new_l)
            this.TV = new_l*this.TV.cost2;
            this.cost = this.F0 + this.TV;
        end
        
        function setBounds(this,new_b)
            this.bounds = new_b;
            this.C = CostRectangle(this.F0.sizein,this.bounds(1),this.bounds(2));
        end
        
        function run(this,x0)
            
            %F=P;
            
            % Reimplementation from :class:`Opti`. For details see [1].
            if ~isempty(x0) % To restart from current state if wanted
                this.xopt = x0;
            end
            if isempty(this.P)
                this.P = zeros(this.D.sizeout);%dual variable
            end
            assert(~isempty(this.xopt),'Missing starting point x0');
            tstart=tic;
            this.OutOp.init();
            this.niter=1;
            this.starting_verb();
            
            if this.Nesterov
                t = 1;
            end
            
            F = this.P;
            lambda = this.TV.cost1;
            
            while (this.niter < this.maxiter)
                
                xold = this.xopt;
                Pnew = F + (1/(this.L*lambda))*this.OpTV*(this.C.applyProx(this.F0.y - lambda*this.OpTV'*(F),0));
                Pnew = Pnew./repmat(max(1,sqrt(sum(Pnew.^2, this.ndims + 1))),[ones(1,this.ndims), this.ndims]);%Project L2 ball
                
                if this.Nesterov
                    tnew = (1 + sqrt(1 + 4*t^2))/2;
                    F = Pnew + (t - 1)/tnew*(Pnew - this.P);
                    t = tnew;
                else
                    F = Pnew;
                end
                this.P = Pnew;
                this.xopt = this.C.applyProx(this.F0.y - lambda*this.OpTV'*(this.P),0);
                
                this.niter = this.niter + 1;
                if this.test_convergence(xold), break; end
				% - Call OutputOpti object
				if (mod(this.niter,this.ItUpOut)==0),this.OutOp.update(this);end
            end
            
            this.xopt = this.C.applyProx(this.F0.y - lambda*this.OpTV'*(this.P),0);
            this.time=toc(tstart);
            this.ending_verb();
        end
    end
end

