classdef OptiNAME < Opti
    %% OptiNAME : TODO ADD DESCRIPTION
    %  Matlab inverse Problems Library
    %
    % -- Description
    % TODO ADD A DESCRIPTION
    %
    % -- Example
    % TODO ADD INSTANTIATION EXAMPLE
    % 
    % -- Properties
    % TODO ADD NEW PROPERTIES
    %
    % -- Methods
    % TODO ADD NEW METHODS
    %
    % -- References
    % TODO ADD REFERENCES IF NEEDED
    %
    % Please refer to the OPTI superclass for general documentation about optimization class
    % See also Opti, OutputOpti
    %
    %     Copyright (C) TODO YEAR NAME EMAIL
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

    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		%% TODO SET HERE NEW PROTECTED SET AND PUBLIC READ PROPERTIES
		%% IF NEEDED.
		%% EXAMPLE THE MINIMIZED FUNC
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
		%% TODO SET HERE NEW FULLY PROTECTED PROPERTIES 
		%% (E.G. INTERNAL VARIABLE USED TO AVOID MULTIPLE COMPUTATION)
    end
    % Full public properties
    properties
		%% TODO SET FULLY PUBLIC PROPERTIES
    end
    
    methods
    	%% Constructor
    	function this=OptiNAME(~)
    		% TODO SET THE INHERITED PROPERTIES
    		this.name='Opti NAME';
    		this.cost=????;
    		this.OutOp=?????;
    		% TODO SET NEW DEFINED PROPERTIES
    	end 
    	%% Run the algorithm
        function run(this,x0) 
			if ~isempty(x0)   % To restart from current state if wanted
				this.xopt=x0;
				% + other initialization that should be avoided for restarting
			end; 
			assert(~isempty(this.xopt),'Missing starting point x0');
			tstart=tic;
			this.OutOp.init();
			this.niter=1;
			this.starting_verb();
			while (this.niter<this.maxiter)
				this.niter=this.niter+1;
				xold=this.xopt;
				% - Algorithm iteration
				
				% TODO IMPLEMENTS THE ALGORITHM ITERATION
				
				% - Convergence test
				if this.test_convergence(xold), break; end
				% - Call OutputOpti object
				if (mod(this.niter,this.ItUpOut)==0),this.OutOp.update(this);end
			end 
			this.time=toc(tstart);
			this.ending_verb();
        end
	end
end
