classdef LinOp
    % LinOP : Linear Operator generic class
    % The LinOp meta class implement generic methods for all linear
    % operators.
    % 
    
    properties (SetAccess = protected,GetAccess = public)
        name = 'none'
        dimsin;
        dimsout;
    end
    
    methods % default methods
        function apply(~,~) 
            error('Apply not implemented');
        end
            function transpose(~,~)
            error('transpose not implemented');
        end
            function y = HtH(self,x) 
            y = self.transpose(self.apply(x));
        end
            function inverse(~,~)
            error('inverse not implemented')
        end        
    end
end

