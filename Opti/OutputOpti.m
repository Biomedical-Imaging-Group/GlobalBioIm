classdef OutputOpti < handle
    %% OutputOpti : Generic class for object OutputOpti objects
    %  Matlab inverse Problems Library
    %
    % At each ItUpOut iterations of an optimization algorithm (see Opti generic class),
    % the update method of an OutputOpti object will be called in order to acheive user 
    % defined computations:
    %   - compute cost / SNR
    %   - store current iterate / cost value 
    %   - plot/display stuffs
    %
    % The present generic class implements a basic update method that:
    %   - display the iteration number
    %   - computes & display the cost (if activated)
    %   - computes & display the SNR if ground truth is provided
    %
    % --Example
    %  VU=OutputOpti(computecost,xtrue,iterVerb)
    % where computecost is a boolean stating if the cost function has to be computed, xtrue is
    % the ground truth and iterVerb corresponds to the interval (in number of iterations) for 
    % which the update will only perform the computation (every ItUpOut iterations) but not
    % print/display anithing. Hence iterVerb must be a multiple of ItUpOut Opti parameter.
    % 
    % Child classes can derive from this one to define different actions to execute during 
    % optimization updates.
    %
    % IMPORTANT: The update method should have an unique imput that is the OPTI object in order to 
    % be generic for all Optimization routines. Hence the update method has acces (in reading mode) 
    % to all the properties of OPTI objects.
    %
    % -- Properties
    % * |name|        - name of the VerbUpdate class
    % * |computecost| - Boolean, if true the cost function will be computed
    % * |xtrue|       - Ground Truth to compute the error with the solution (if provided)
    % * |evolcost|    - array to save the evolution of the cost function
    % * |evolsnr|     - array to save the evolution of the SNR
    % * |iterVerb|    - message will be displayed every iterVerb iterations (must be a multiple
    %                   of the ItUpOut Opti parameter)
    %
    % -- Methods
    % * |update|  - execute the defined actions
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
        name = 'OutputOpti'% name of the optimization algorithm
		computecost=false; % Boolean, if true the cost function will be computed
		xtrue;             % Ground Thruth
	    evolcost;          % array saving the evolution of the cost function
		evolsnr;           % array saving the evolution of the error with the groud truth
		evolxopt;          % cell saving the optimization variable xopt
		iternum;           % array saving the iteration number corresponding to evolcost, evolxopt and evolerr entries
		iterVerb=0;        % message will be displayed every iterVerb iterations
    end
    properties (SetAccess = protected,GetAccess = public)
    	normXtrue;         % Norm of the true signal (to compute snr)
    	isgt=false;        % Boolean true if Ground Truth is provided
   		count;             % internal counter
    end
    
    methods
    	%% Constructor
        function this=OutputOpti(computecost,xtrue,iterVerb) 
        	if nargin==1
        		this.computecost=computecost;
            elseif nargin==2
            	this.computecost=computecost;
            	this.xtrue=xtrue;
            elseif nargin==3
                this.computecost=computecost;
            	this.xtrue=xtrue;
            	this.iterVerb=iterVerb;
			end
			if ~isempty(this.xtrue)
            	this.isgt=true;
            	this.xtrue=xtrue;
            	this.normXtrue=norm(this.xtrue(:));
            end
        end
        %% Initialization
        function init(this)
        	this.count=1;
        	this.evolcost=[];
        	this.evolsnr=[];
			this.iternum = [];
			this.evolxopt = {};
        end
        %% Update method
        function update(this,opti)
        	str=sprintf('Iter: %5i',opti.niter);
        	if this.computecost
        		cc=opti.cost.apply(opti.xopt);
        		str=sprintf('%s | Cost: %4.4e',str,cc);
        		this.evolcost(this.count)=cc;
        	end
        	if this.isgt
        		snr=20*log10(this.normXtrue/norm(this.xtrue(:)-opti.xopt(:)));
        		str=sprintf('%s | SNR: %4.4e dB',str,snr);
        		this.evolsnr(this.count)=snr;
        	end
        	this.evolxopt{this.count}=opti.xopt;
        	this.iternum(this.count)=opti.niter;
        	this.count=this.count+1;
        	if (mod(opti.niter,this.iterVerb)==0) || opti.niter==1,
        		disp(str);
        	end
        end
    end
end
