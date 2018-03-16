classdef LinOpShape <  LinOp
    % LinOpShape: reshaping operator
    %
    % Reshape an array of size sizein in a array of size sizeout
    %
    % :param sizein: input size
    % :param sizeout: output size
    %
    % **Example** R=LinOpShape(sizein,sizeout)
    %
    % See also :class:`LinOp`, :class:`Map`
    
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
    
    %% Constructor
    methods
        function this = LinOpShape(sizein, sizeout)
            this.name ='LinOpShape';
            this.isInvertible=true;
			
			this.norm = 1;
            
            assert(issize(sizein),'The input size sizein should be a conformable  to a size ');
            this.sizein = sizein;
            assert(issize(sizeout),'The input size sizeout should be a conformable  to a size ');
            this.sizeout = sizeout;
            
            assert(prod(sizeout)==prod(sizein),'The inputs sizein and sizeout should be a conformable');
            
		end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.     
            y = reshape(x, this.sizeout);
        end	
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.     
            y = reshape(x, this.sizein);
        end		
        function y = applyHHt_(~,x)
            % Reimplemented from parent class :class:`LinOp`.     
            y=x;
        end		
        function y = applyHtH_(~,x)
            % Reimplemented from parent class :class:`LinOp`.     
            y=x;
		end   		
        function y = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.     
             y = reshape(x, this.sizein);
        end		
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.     
            y = reshape(x, this.sizeout);
        end
        
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpIdentity(this.sizeout);
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpIdentity(this.sizein);
        end
    end
end
