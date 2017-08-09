classdef LinOpIdentity <  LinOp
    % Identity linear operator
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{x}$$
    %
    % :param sz: size of \\(\mathrm{x}\\) on which the :class:`LinOpIdentity` applies.
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %%    Copyright (C) 2015 
    %     F. Soulez  ferreol.soulez@epfl.ch
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
	  function this = LinOpIdentity(sz)
		this.name='LinOp Identity';
		this.isComplexIn=true;
		this.isComplexOut=true;
        this.isDifferentiable=true;
        this.isInvertible=true;
		this.norm=1;
		if nargin>0
		  this.sizeout=sz;
		  this.sizein=sz;
		else
		  error('A size SZ should be given');
		end
	  end
	end
	
	%% Core Methods containing implementations (Protected)
    methods (Access = protected)      
        function y = apply_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.       	
            y =x;
        end        
        function y = applyAdjoint_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.      	
            y =x;
        end        
        function y = applyHtH_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.       	
            y =x;
        end       
        function y = applyHHt_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.       	
            y =x;
        end       
        function y = applyInverse_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.        	
            y =x;
        end        
        function y = applyAdjointInverse_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.        	
            y =x;
        end
        function M = makeComposition_(this,G)
            % Reimplemented from parent class :class:`LinOp`.
            % Returns \\(\\mathrm{G}\\).
            M = G;
        end
    end
end

