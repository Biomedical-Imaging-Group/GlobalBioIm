classdef VerbUpdate < handle
    %% VerbUpdate : Generic class for object VerbUpdates objects 
    %  Matlab Inverse Problems Library
    %
    % At each verb iterations of an optimization algorithm (see Opti generic class),
    % a function can be executed in order to:
    %   - compute cost / error with ground thruth
    %   - store current iterate / cost value 
    %   - plot/display stuffs
    %
    % The present generic class implements a basic exec method that:
    %   - display the iteration number
    %   - computes & display the cost (if activated)
    %   - computes & display the error to the groung thruth is provided
    %
    % --Example
    %  VU=VerbUpdate(computecost,xtrue)
    % 
    % Child classes can derive from this one to define different actions to execute during 
    % optimization updates.
    %
    % IMPORTANT: The exec method should have an unique imput that is the OPTI object in order to 
    % be generic for all Optimization routines. Hence the exec method has acces (in reading mode) 
    % to all the properties of OPTI objects.
    %
    % -- Properties
    % * |name|        - name of the VerbUpdate class
    % * |computecost| - Boolean, if true the cost function will be computed
    % * |xtrue|       - Ground Thruth to compute the error with the solution (if provided)
    % * |evolcost|    - array to save the evolution of the cost function
    % * |evolerr|     - array to save the evolution of the error with ground thruth
    %
    % -- Methods
    % * |exec|    - execute the defined actions
    % * |init|    - initialization (called at the starting of the opti algorithm to initialize
    %               internal variables (counter, array for saving...)
    %
    % See also Opti
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
    
    properties (SetAccess = protected,GetAccess = public)
        name = 'none'      % name of the optimization algorithm
        isgt=false;        % Boolean true if Ground Thruth is provided
		computecost=false; % Boolean, if true the cost function will be computed
		xtrue;             % Ground Thruth
	    count;             % internal counter
	    evolcost;          % array saving the evolution of the cost function
		evolerr;           % array saving the evolution of the error with the groud thruth
		evolxopt;          % cell saving the optimization variable xopt
		iternum;           % array saving the iteration number corresponding to evolcost, evolxopt and evolerr entries
    end
    
    methods
    	%% Constructor
        function this=VerbUpdate(computecost,xtrue) 
        	if nargin==1
        		this.computecost=computecost;
            elseif nargin==2
            	this.computecost=computecost;
            	this.isgt=true;
            	this.xtrue=xtrue;
			end
        end
        %% Initialization
        function init(this)
        	this.count=1;
        	this.evolcost=[];
        	this.evolerr=[];
        end
        %% Exec method
        function exec(this,opti)
        	str=sprintf('Iter: %5i',opti.niter);
        	if this.computecost
        		cc=opti.cost.eval(opti.xopt);
        		str=sprintf('%s | Cost: %4.4e',str,cc);
        		this.evolcost(this.count)=cc;
        	end
        	if this.isgt
        		err=norm(this.xtrue(:)-opti.xopt(:));
        		str=sprintf('%s | GT-Err: %4.4e',str,err);
        		this.evolerr(this.count)=err;
        	end
        	this.evolxopt{this.count}=opti.xopt;
        	this.iternum(this.count)=opti.niter;
        	this.count=this.count+1;
        	disp(str);
        end
    end
end
