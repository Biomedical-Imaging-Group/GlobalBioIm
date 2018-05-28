classdef OutputOptiConjGrad < OutputOpti
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
