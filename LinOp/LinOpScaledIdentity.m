classdef LinOpScaledIdentity <  LinOp
    % Identity linear operator
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\nu\\mathrm{x}$$
    % where \\(\\nu \\in \\mathbb{R}\\).
    %
    % :param sz: size of \\(\mathrm{x}\\) on which the :class:`LinOpIdentity` applies.
    % :param nu: scaling parameter (default 1)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % See also :class:`LinOp`, :class:`Map`
    
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
    
    properties
        nu=1;
    end
    
	%% Constructor
	methods
	  function this = LinOpScaledIdentity(sz,nu)
		this.name='LinOp ScaledIdentity';
        if nargin==2
            this.nu=nu;
        end
		this.isComplexIn=true;
		this.isComplexOut=true;
        this.isDifferentiable=true;
        if this.nu~=0
            this.isInvertible=true;
        end
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
            y =this.nu*x;
        end        
        function x = applyAdjoint_(this,y)
        	% Reimplemented from parent class :class:`LinOp`.      	
            x =this.nu*y;
        end        
        function y = applyHtH_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.       	
            y =this.nu^2*x;
        end       
        function x = applyHHt_(this,y)
        	% Reimplemented from parent class :class:`LinOp`.       	
            x =this.nu^2*y;
        end       
        function x = applyInverse_(this,y)
        	% Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                x =y/this.nu;
            end
        end        
        function y = applyAdjointInverse_(this,x)
        	% Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                y =x/this.nu;
            end
        end
    end
end