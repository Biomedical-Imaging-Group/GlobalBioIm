classdef Opti < matlab.mixin.SetGet
    % Abstract class for optimization algorithms to minimize :class:`Cost` objects
    %
    % :param name: name of the algorithm
    % :param cost: minimized :class:`Cost`
    % :param maxiter: maximal number of iterations (default 50)
    % :param verbose: bollean (default true) to activate verbose mode
    % :param OutOp: :class:`OutputOpti` object
    % :param ItUpOut: number of iterations between two calls to the update method of the  :class:`OutputOpti` object :attr:`OutOp` (default 0)
    % :param CvOp: :class:`TestCvg` object
    % :param time: execution time of the algorithm
    % :param niter: iteration counter
    % :param xopt: optimization variable
    %
    % See also :class:`OutputOpti` :class:`Cost`
    
    %%    Copyright (C) 2017
    %     E. Soubies emmanuel.soubies@epfl.ch
    %     F. Soulez ferreol.soulez@univ-lyon1.fr
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
    
    properties (Constant)
        OPTI_NEXT_IT = 0;   % new iteration
        OPTI_REDO_IT = 1; %
        OPTI_STOP    = 2;    % stop iteration
    end
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        name = 'none'        % name of the optimization algorithm
        cost;                % minimized cost
        time;                % running time of the algorithm (last run)
        niter;               % iteration counter
        xopt=[];             % optimization variable
        xold=[];
        needxold=false;
    end
    % Full public properties
    properties
        endingMessage;       % Ending message
        verbose=true;        % if true display information (starting and ending message
        OutOp=OutputOpti(false,1);  % OutputOpti object
        CvOp=TestCvg();      % OutputOpti object
        maxiter=50;     % maximal number of iterates
        ItUpOut=0;      % period (in number of iterations) of calling the OutputOpti object
    end
    
    %% Methods
    methods
        function run(this,x0)
            % Run the algorithm.
            %
            % :param x0: initial point in \\(\\in X\\), if no argument restarts from the current value :attr:`xopt`.
            %
            % **note**: this method does not return anything, the result being stored in public attribute :attr:`xopt`.
            if(nargin==1)
                assert(~isempty(this.xopt),'Missing starting point x0');
            else
                this.initialize(x0);
            end            
            tstart=tic;
            this.OutOp.init();
            this.niter=0;
            this.starting_verb();
            this.endingMessage= ['Maximum number of iterations reached: ', num2str(this.maxiter)];           
            while (this.niter<this.maxiter)
                % - Update parameters
                this.updateParams();
                % - Algorithm iteration   
                flag=this.doIteration();

                if(flag == this.OPTI_NEXT_IT )
                    this.niter=this.niter+1;
                    
                    % - Convergence test
                    if this.niter>1 && this.CvOp.testConvergence(this), break; end
                    
                    % - Call OutputOpti object
                    if this.ItUpOut>0 && (((mod(this.niter,this.ItUpOut)==0)|| (this.niter==1)))
                        this.OutOp.update(this);
                    end
                    if this.needxold
                        this.xold=this.xopt;
                    end
                elseif  flag==this.OPTI_STOP,
                    break;
                end
            end
            this.time=toc(tstart);
            this.ending_verb();
        end
        function initialize(this,x0)
            % Implements initialization of the algorithm
            %
            % :param x0: initial point
            
            if ~isempty(x0) % To restart from current state if wanted
                this.xopt=x0;
                if this.CvOp.needxold
                    this.needxold=true;
                end
                if this.needxold
                    this.xold=x0;
                end
            end
        end
        function flag=doIteration(this)
            % Implements algorithm iteration
            %
            % :return: flag with values
            %
            %    - OPTI_NEXT_IT (= 0) to go to the next iteration
            %    - OPTI_REDO_IT (= 1) to redo the iteration
            %    - OPTI_STOP    (= 2) to stop algorithm
            
            error(['In ',this.name,': doIteration method is not implemented']);
        end
        function updateParams(this)
            % Updates the parameters of the algorithm at each iteration
            % (default: no update). This method can be overloaded to makes
            % some parameters varying during iterations (e.g. descent step,
            % lagrangian parameters...)
        end
        function starting_verb(this)
            % Generic method to display a starting message in verbose mode.
            
            if this.verbose
                fprintf('---> Start %s ... \n',this.name);
                
            end
        end
        function ending_verb(this)
            % Generic method to display a ending message in verbose mode.
            
            if this.verbose
                disp(this.endingMessage);
                fprintf('... Optimization finished \nElapsed time (s): %4.2d (%i iterations). \n',this.time, this.niter);
            end
        end
    end
end
