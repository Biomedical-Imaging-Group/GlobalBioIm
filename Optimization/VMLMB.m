classdef VMLMB
    properties (Constant)
        OP_TASK_START  = 0; % first entry, start search
        OP_TASK_FG     = 1; % computation of F and G requested
        OP_TASK_NEWX   = 2; % new improved solution available for inspection
        OP_TASK_CONV   = 3; % search has converged
        OP_TASK_WARN   = 4; % search aborted with warning
        OP_TASK_ERROR  = 5; % search aborted with error
    end
    properties
        m=3;
        fatol=0.0;
        frtol=1e-12;
        sftol=0.001;
        sgtol=0.9;
        sxtol=0.1;
        epsilon=0.01;
        costheta=0.4;
        nbitermax;
        nbevalmax;
    end
    properties (SetAccess = protected,GetAccess = public)
        nparam;
        nbeval=0;
        task;
        bounds;
        xmin;
        xmax;
    end
    methods
        function this = VMLMB(nparam,xmin,xmax)
            this.nparam =nparam;
            if(nargin>1)
            this.xmin = xmin;
            this.xmax = xmax;
            end
        end
        
end