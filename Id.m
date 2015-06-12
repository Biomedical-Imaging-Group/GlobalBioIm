classdef Id <  LinOp
    % LinOP : Linear Operator generic class
    % The LinOp meta class implement generic methods for all linear
    % operators.
    % 
    % 
    
    methods % default methods
        function y = apply(~,x) 
            y =x;
        end
            function transpose(~,x)
            y =x;
        end
            function y = HtH(~,x) 
            y =x;
        end
            function inverse(~,x)
            y =x;
        end        
%         function dimsin = get.dimsin(self)
%         end
%         function dimsout = get.dimsout(self)
%         end
%         function name = get.name(self)
%         end
    end
end

