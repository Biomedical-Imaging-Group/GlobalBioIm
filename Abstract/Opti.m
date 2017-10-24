classdef Opti < matlab.mixin.SetGet
    % Abstract class for optimization algorithms to minimize :class:`Cost` objects
    %
    % :param name: name of the algorithm
    % :param cost: minimized :class:`Cost`
    % :param maxiter: maximal number of iterations (default 50)
    % :param xtol: tolerance on the relative difference between two iterates (default 1e-5)
    % :param OutOp: :class:`OutputOpti` object
    % :param ItUpOut: number of iterations between two calls to the update method of the  :class:`OutputOpti` object :attr:`OutOp` (default 0)
	% :param time: execution time of the algorithm
	% :param niter: iteration counter
	% :param xopt: optimization variable
    %
    % See also :class:`OutputOpti` :class:`Cost`

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

    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
        name = 'none'        % name of the optimization algorithm
    	cost;                % minimized cost
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
        function run(this,x0) 
        	% Run the algorithm. 
        	%
        	% :param x0: initial point in \\(\\in X\\), if x0=[] restarts from the current value :attr:`xopt`.
        	%
			% **note**: this method does not return anything, the result being stored in public attribute :attr:`xopt`.

            error(['In ',this.name,': run method is not implemented']);
        end
        function starting_verb(this)  
        	% Generic method to display a starting message in verbose mode.
        	   	
            if this.ItUpOut~=0
                if this.OutOp.iterVerb~=0
                    fprintf('---> Start %s ... \n',this.name);
                end
                this.OutOp.update(this);
            end
        end
        function ending_verb(this)
        	% Generic method to display a ending message in verbose mode.
        	
        	if this.ItUpOut~=0 && this.OutOp.iterVerb~=0 
				fprintf('... Optimization finished \nElapsed time (s): %4.2d (%i iterations). \n',this.time, this.niter);
        	end
        end
        function stop=test_convergence(this,xold)
        	% Tests algorithm convergence from the relative difference between two successive iterates 
        	%
        	% :param xold: iterate \\( \\mathrm{x}^{k-1}\\).
        	% :returns stop: boolean true if
        	% $$ \\frac{\\| \\mathrm{x}^{k} - \\mathrm{x}^{k-1}\\|}{\\|\\mathrm{x}^{k-1}\\|} < x_{tol}.$$
        	
        	r=this.xopt-xold;
        	xdiff=norm(r(:))/(norm(xold(:))+eps);
        	stop=xdiff<this.xtol;
        end
    end
end
