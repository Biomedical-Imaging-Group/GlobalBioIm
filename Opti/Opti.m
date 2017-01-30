classdef Opti < handle
    %% Opti : Optimization algorithm generic class
    %  Matlab Inverse Problems Library
    %  The Opti meta class implements generic methods for all optimization algorithms
    %
    % -- Properties
    % * |name|      - name of the optimization algorithm  
    % * |maxiter|   - maximal number of iterations (default 50) public property
    % * |xtol|      - stopping criteria tolerance on the relative difference between two 
    %                 iterates (see TestConvergence function, default 1e-5) public property
    % * |ItUpOut|   - every ItUpOut iterations the update method of a OutputOpti object is called 
    %                 (see Class OutputOpti, default 0) public property
    % * |OutOp|     - OutputOpti object
    % * |cost|      - whole minimized functional (Func)
    % * |time|      - execution time of the algorithm (last run)
    % * |niter|     - iteration counter
    % * |xopt|      - optimization variable
    %
    % Note: in each derived class the minimized functional(s) are defined as properties and combined to 
    %       form the cost functional (that can be used by the update method of OutputOpti object).
    %       For example let us consider the OptiChambollePock class that allows to minimize F(Kx)+G(x). In this 
    %       situation F, G and the LinOp K have to be given separately to the Opti class (they will define
    %       new properties of the derivate class). Moreover, we will define cost = F(K) + G in order to have a generic
    %       way to evaluate the cost for all the algorithms (which can be used by the method update of OutputOpti)
    %
    % -- Methods
    % * |run(x0)|          - run the algorithm from the initial point x0. If x0=[], should restart from current state. 
    % * |starting_verb|    - generic method to display a starting message in verb mode
    % * |ending_verb|      - idem for ending message
    % * |test_convergence| - generic function to test convergence based on the relative difference between two successive iterates
    %
    % See also OutputOpti
    %
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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
        name = 'none'        % name of the optimization algorithm
    	cost;                % minimized functional
    	time;                % running time of the algorithm (last run)
    	niter;               % iteration counter
    	xopt=[];             % optimization variable
    	OutOp=OutputOpti();  % OutputOpti object
    end
    % Full public properties
    properties
    	maxiter=50;     % maximal number of iterates
    	xtol=1e-5;      % stopping criteria tolerance on the relative difference btw two successive iterates
    	ItUpOut=0;      % period (in number of iterations) of calling the OutputOpti object
    end
    
    methods
    	%% Run the algorithm
        function run(~,~) 
            error(['In ',this.name,': run method is not implemented']);
        end
        %% Display starting message
        function starting_verb(this)     	
        	if this.ItUpOut~=0
        		fprintf('---> Start %s ... \n',this.name);
				this.OutOp.update(this);
			end
        end
        %% Display ending message
        function ending_verb(this)
        	if this.ItUpOut~=0
				fprintf('... Optimization finished \nElapsed time (s): %4.2d (%i iterations). \n',this.time, this.niter);
        	end
        end
        %% Test Convergence with the relative difference btw two successive iterates
        function stop=test_convergence(this,xold)
        	r=this.xopt-xold;
        	xdiff=norm(r(:))/(norm(xold(:))+eps);
        	stop=xdiff<this.xtol;
        end
    end
end
