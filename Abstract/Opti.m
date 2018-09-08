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
 
   
    %% Properties
    % - Public 
    properties (SetObservable, AbortSet)
        verbose=true;        % if true display information (starting and ending message
        OutOp=OutputOpti(false,[],1);  % OutputOpti object
        CvOp=TestCvg();      % OutputOpti object
        maxiter=50;          % maximal number of iterates
        ItUpOut=0;           % period (in number of iterations) of calling the OutputOpti object
        cost;                % minimized cost
        endingMessage;       % Ending message
    end
    % - Readable
    properties (SetAccess = protected,GetAccess = public)
        name = 'none'        % name of the optimization algorithm
        time;                % running time of the algorithm (last run)
        niter;               % iteration counter
        xopt=[];             % optimization variable
        xold;
    end
    % - Protected
    properties (Access = protected)       
        listenerList=cell(0);
    end
    % - Hidden
    properties (Hidden)
        cleanup
    end
    % - Constant
    properties (Hidden,Constant)
        OPTI_NEXT_IT = 0;   % new iteration
        OPTI_REDO_IT = 1;   %
        OPTI_STOP    = 2;   % stop iteration
    end
    
    %% Events
    events
        modified % to propagate the modification of a property into compositions
    end
    
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Implements initializing actions for a given property (prop)
        end
        function initObject(this,hrch)
            % Run the updateProp method with parameter 'all' and
            % initialize the listeners if this method is executed from the
            % hierarchy level hrch corresponding to class(this).
            if strcmp(class(this),hrch)
                this.updateProp('all');
                prop=properties(this);
                for i=1:length(prop)
                    % If public property, add a PostSet listener
                    pp=findprop(this,prop{i});
                    if strcmp(pp.SetAccess,'public')
                        this.listenerList{end+1}=addlistener(this,prop{i},'PreSet',@this.handlePreSetEvents);
                        this.listenerList{end+1}=addlistener(this,prop{i},'PostSet',@this.handlePostSetEvents);
                        % If property is a Map, add a modified listener
                        if isa(this.(prop{i}),'Map')
                            this.listenerList{end+1}=addlistener(this.(prop{i}),'modified',@this.handleModifiedEvents);
                        elseif isa(this.(prop{i}),'cell') && isa(this.(prop{i}){1},'Map')
                            for j=1:length(this.(prop{i}))
                                this.listenerList{end+1}=addlistener(this.(prop{i}){j},'modified',@this.handleModifiedEvents);
                            end
                        end
                    end
                end
                this.cleanup = onCleanup(@()delete(this));
            end
        end
        function handlePostSetEvents(this,src,evnt)
            % This function is called when an observable property is
            % set to another value in order to update properly some
            % (precomputed) properties (re-initialize)
            
            %disp(['     In ',this.name,' property ',src.Name,' has been updated']);
            this.updateProp(src.Name);
            % Because the listeners for modified events are lost when a new
            % object is set to the property -> define a new one
            if strcmp(evnt.EventName,'PostSet')
                if isa(this.(src.Name),'Map')
                    this.listenerList{end+1}=addlistener(this.(src.Name),'modified',@this.handleModifiedEvents);
                elseif isa(this.(src.Name),'cell') && isa(this.(src.Name){1},'Map')
                    for i=1:length(this.(src.Name))
                        if ~any(cellfun(@(x) isequal(this.(src.Name){i},x.Source{1}),this.listenerList))
                            this.listenerList{end+1}=addlistener(this.(src.Name){i},'modified',@this.handleModifiedEvents);
                        end
                    end
                end
            end
            % To propagate the modification within Map properties
            notify(this,'modified');
        end
        function handlePreSetEvents(this,src,~)
            % This function is called before an observable property is set
            if isa(this.(src.Name),'Map')
                this.(src.Name).clearListenerList();
            elseif isa(this.(src.Name),'cell') && isa(this.(src.Name){1},'Map')
                for i=1:length(this.(src.Name))
                    this.(src.Name){i}.clearListenerList();
                end
            end
        end
        function handleModifiedEvents(this,src,evnt) % Necessary for properties which are objects of the Library (i.e. Maps)
            % This function is called when an observable property has been
            % modified in order to update properly some (precomputed)
            % properties (re-initialize)
            prop=properties(this);
            if isvalid(this) % to avoid error with already deleted objects
                idx=cell2mat(cellfun(@(x) (isa(this.(x),'Map')  && this.(x)==src) || ...
                    (isa(this.(x),'cell') && any(cell2mat(cellfun(@(y) isa(y,'Map') && y==src,this.(x),'UniformOutput',false))))...
                    ,prop,'UniformOutput',false));
                if any(idx)
                    sourc.Name=prop{idx};
                    %disp(['     In ',this.name,' property ',sourc.Name,' has been modified']);
                    handlePostSetEvents(this,sourc,evnt);
                end
            end
        end
    end
    
    %% Methods for optimization
    % - run
    % - initialize
    % - doIteration
    % - updateParams
    % - starting_verb
    % - ending_verb
    methods
        function run(this,x0)
            % Run the algorithm.
            %
            % :param x0: initial point in \\(\\in X\\), if no argument restarts from the current value :attr:`xopt`.
            %
            % **note**: this method does not return anything, the result being stored in public attribute :attr:`xopt`.
            if(nargin==1)
                assert(~isempty(this.xopt),'Missing starting point x0');
                x0 = this.xopt;
            end
            this.initialize(x0);
            tstart=tic;
            this.OutOp.init();
            this.niter=0;
            this.starting_verb();
            this.endingMessage= ['Maximum number of iterations reached: ', num2str(this.maxiter)];
            this.xold= x0;
            while (this.niter<this.maxiter)
                % - Update parameters
                this.updateParams();
                % - Algorithm iteration   
                flag=this.doIteration();

                if(flag == this.OPTI_NEXT_IT)
                    this.niter=this.niter+1;
                    
                    % - Convergence test
                    if this.CvOp.testConvergence(this), break; end
                    
                    % - Call OutputOpti object
                    if this.ItUpOut>0 && (((mod(this.niter,this.ItUpOut)==0)|| (this.niter==1)))
                        this.OutOp.update(this);
                    end
                    this.xold=this.xopt;
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
    
    %% Utility methods
    % - clearListenerList
    methods
     function clearListenerList(this,obj)
            % Clear properly the listener list and call recursively on all
            % public properties.
            % If the given object is this, then do nothing
            if nargin==1, obj={}; end;
            if ~any(cellfun(@(x) isequal(this,x),obj))
                prop=properties(this);
                for i=1:length(prop)
                    % If public property, add a PostSet listener
                    pp=findprop(this,prop{i});
                    if ~isempty(pp) && strcmp(pp.SetAccess,'public')
                        if isa(this.(prop{i}),'Map')
                            this.(prop{i}).clearListenerList(obj);
                        elseif isa(this.(prop{i}),'cell') && isa(this.(prop{i}){1},'Map')
                            for j=1:length(this.(prop{i}))
                                this.(prop{i}){j}.clearListenerList(obj);
                            end
                        end
                    end
                end
                if isvalid(this) % to avoid error with already deleted objects
                    for ii = 1:numel(this.listenerList)
                        if isa(this.listenerList{ii}.Source{1},'Map')
                            delete(this.listenerList{ii});
                        end
                    end
                    this.listenerList=cell(0);
                end
            end
        end
    end
end
