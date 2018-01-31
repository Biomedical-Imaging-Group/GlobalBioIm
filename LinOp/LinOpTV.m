classdef LinOpTV < LinOp
    properties
        D;
        ndms;
        bc;
        allElements;
    end
    methods
        function this = LinOpTV(sz,bc)
            if nargin <=1 || isempty(bc), bc = 'circular';end
            this.name = 'LinOpTV';
            this.isDifferentiable = false;
            
            this.bc = bc;
            this.ndms = numel(sz);
            this.allElements = repmat({':'}, 1, this.ndms);
            this.sizein = sz;
            this.sizeout = [sz,this.ndms];
            
            for n = 1:this.ndms
                this.D{n} = LinOpGrad(sz,n,this.bc);
            end
        end
    end
        %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
    % - applyHtH_(this,x)
    % - makeHtH_(this)
    methods (Access = protected)
        function y = apply_(this,x)
            y = zeros(this.sizeout);
            for n = 1:this.ndms
                y(this.allElements{:},n) = this.D{n}*x;
            end
        end
        
        function y = applyAdjoint_(this,x)
            y = zeros(this.sizein);
            for n = 1:this.ndms
                y = y + (this.D{n}')*x(this.allElements{:},n);
            end
            
        end
    end
end