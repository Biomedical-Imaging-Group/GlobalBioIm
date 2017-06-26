
classdef LinOpIdentity <  LinOp
    % Identity operator
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{x}$$
    %
    % :param sz: size of \\(\mathrm{x}\\) on which the :class:`LinOpIdentity` applies.
    %
    % See also :class:`LinOp`
    
    %     Copyright (C) 2015 F. Soulez  ferreol.soulez@epfl.ch
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
    
    methods
        function this = LinOpIdentity(sz)
            this.name='LinOp Identity';
            this.iscomplex=false;
            this.isinvertible=true;
            this.norm=1;
            if nargin>0
                this.sizeout=sz;
                this.sizein=sz;
            else
                error('a size SZ should be given');
            end
        end
        
        function y = apply(this,x)
        	% Reimplemented from parent class :class:`LinOp`.
        	
            y =x;
        end
        
        function y = adjoint(this,x)
        	% Reimplemented from parent class :class:`LinOp`.
        	
            y =x;
        end
        
        function y = HtH(this,x)
        	% Reimplemented from parent class :class:`LinOp`.
        	
            y =x;
        end
        
        function y = HHt(this,x)
        	% Reimplemented from parent class :class:`LinOp`.
        	
            y =x;
        end
        
        function y = inverse(this,x)
        	% Reimplemented from parent class :class:`LinOp`.
        	
            y =x;
        end
        
        function y = adjointInverse(this,x)
        	% Reimplemented from parent class :class:`LinOp`.
        	
            y =x;
        end
    end
end

