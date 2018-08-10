classdef (Abstract) Map < handle
    % Abstract class for Maps which maps elements from \\(\\mathrm{X}\\) to
    % \\(\\mathrm{Y}\\)
    % $$ \\mathrm{H}: \\mathrm{X}\\rightarrow \\mathrm{Y}.$$
    % where \\(\\mathrm{X}\\) and \\(\\mathrm{Y}\\) are either
    % \\(\\mathbb{R}^N\\) or \\(\\mathbb{C}^N\\).
    %
    % :param name: name of the linear operator \\(\\mathbf{H}\\)
    % :param sizein:  dimension of the left hand side vector space \\(\\mathrm{X}\\)
    % :param sizeout:  dimension of the right hand side vector space \\(\\mathrm{Y}\\)
    % :param norm: norm of the operator \\(\\|\\mathrm{H}\\|\\) (if known, otherwise -1)
    % :param isInvertible:  true if the method :meth:`applyInverse_` is implemented
    % :param isDifferentiable:  true if the method :meth:`applyJacobianT_` is implemented
    % :param memoizeOpts: structure of boolean (one field per method, see details below).
    % :param doPrecomputation: boolean true to allow doing precomputations to save time (will generally require more memory).
    %
    % **Note on the memoize option** This option allows to store the result
    % of a method such that if an identical call to this method is done,
    % calculations are avoided. Example: memoizeOpts.apply=true will store
    % the result of H*x.
    
    %%    Copyright (C) 2017
    %     M. McCann michael.mccann@epfl.ch &
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
    
    %% Properties
    % - Public 
    properties (SetObservable, AbortSet)
        memoizeOpts = struct('apply', false, ...
            'applyJacobianT', false, ...
            'applyInverse', false);
        doPrecomputation = false;
    end   
    % - Readable 
    properties (SetAccess = protected)
        name = 'none'             % name of the linear operator
        isInvertible = false;     % true if H.applyInverse(  ) will work
        isDifferentiable = false; % true if H.applyJacobianT(   ) will work
        sizein;                   % dimension of the right hand side vector space
        sizeout;                  % dimension of the left hand side vector space
        norm=-1;                  % norm of the operator
    end
    % - Private
    properties (SetAccess = protected,GetAccess = protected)
        memoCache = struct('apply', struct('in', [], 'out', []),...
            'applyJacobianT', struct('in', [], 'out', []), ...
            'applyInverse', struct('in', [], 'out', []));
        precomputeCache = struct();
        listenerList=cell(0);
    end
    % - Hidden
    properties (Hidden)
        cleanup
    end
    
    %% Events
    events
      modified % to propagate the modification of a property into compositions
    end
        
    %% Initialize and methods related to listeners (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Implements initializing actions for a given property (prop)
            if strcmp(prop,'memoizeOpts') 
                fnames = fieldnames(this.memoizeOpts); % clean the memoize cache
                for i = 1:length(fnames)
                    if ~this.memoizeOpts.(fnames{i})
                        this.memoCache.(fnames{i})=struct('in', [], 'out', []);
                    end
                end
            end
            if strcmp(prop,'doPrecomputation')
                this.precomputeCache = struct();  % clean the precompute cache
            end
            if strcmp(prop,'all') || (~strcmp(prop,'memoizeOpts') && ~strcmp(prop,'doPrecomputation'))
                this.clearCaches();
            end
        end
        function initialize(this,hrch)
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
    
    %% Interface Methods (cannot be overloaded in derived classes: Sealed)
    % - apply(this,x)
    % - applyJacobianT(this, y, v)
    % - applyInverse(this,y)
    % - makeComposition(this, G)
    % - plus(this,G)
    % - minus(this,G)
    % - mpower(this,p)
    % - mtimes(this,G)
    % - times(this,G)
    % - size(this, [dim])
    methods (Sealed)
        function x = apply(this, x)
            % Computes \\(\\mathrm{y}=\\mathrm{H}(\\mathrm{x})\\) for the given
            % \\(\\mathrm{x} \\in \\mathrm{X}\\).
            %
            % Calls the method :meth:`apply_`
            
            if ~checkSize(x, this.sizein) % check input size
                error('Input to apply was size [%s], didn''t match  %s sizein: [%s].',...
                    num2str(size(x)),class(this), num2str(this.sizein));
            end
            % memoize
            if this.memoizeOpts.apply
                x = this.memoize('apply', @this.apply_, x);
            else
                x= this.apply_(x);
            end
            % check output size
            if ~checkSize(x, this.sizeout)
                warning('Output of apply was size [%s], didn''t match %s sizeout: [%s].',...
                    num2str(size(x)),class(this), num2str(this.sizeout));
            end
        end
        function x = applyJacobianT(this, x, v)
            % Compute \\(\\mathrm{x}=[\\mathrm{J}_{\\mathrm{H}}(\\mathrm{v})]^{\\star}\\mathrm{y}\\)
            % where
            %
            % - \\([\\mathrm{J}_{\\mathrm{H}}(\\mathrm{v})]\\) is the Jacobian matrix of
            %   the Map \\(\\mathrm{H}\\) computed at \\(\\mathrm{v} \\in \\mathrm{X} \\)
            % - \\(\\mathrm{y} \\in \\mathrm{Y} \\)
            %
            % Calls the method :meth:`applyJacobianT_`
            
            if ~checkSize(x, this.sizeout) % check input size
                error('Input of y applyJacobianT was size [%s], didn''t match %s sizeout: [%s].',...
                    num2str(size(x)),class(this), num2str(this.sizeout));
            end
            if ~checkSize(v, this.sizein)
                error('Input to v applyJacobianT was size [%s], didn''t match %s sizein: [%s].',...
                    num2str(size(v)),class(this), num2str(this.sizein));
            end
            % memoize
            if this.memoizeOpts.applyJacobianT
                x = this.memoize('applyJacobianT', @this.applyJacobianT_, {x, v});
            else
                x= this.applyJacobianT_(x,v);
            end
            % check output size
            if ~checkSize(x, this.sizein)
                warning('Output of applyJacobianT was size [%s], didn''t match %s sizein: [%s].',...
                    num2str(size(x)),class(this), num2str(this.sizein));
            end
        end
        function x = applyInverse(this, x)
            % Computes \\(\\mathrm{x} = \\mathrm{H}^{-1} \\mathrm{y}\\) for the given
            % \\(\\mathrm{y} \\in \\mathrm{Y}\\). (if applicable)
            %
            % Calls the method :meth:`applyInverse_`
            if ~checkSize(x, this.sizeout) % check input size
                error('Input to applyInverse was size [%s], didn''t match %s sizeout: [%s].',...
                    num2str(size(x)),class(this), num2str(this.sizeout));
            end
            % memoize
            if this.memoizeOpts.applyInverse
                x = this.memoize('applyInverse', @this.applyInverse_,x);
            else
                x= this.applyInverse_(x);
            end
            % check output size
            if ~checkSize(x, this.sizein)
                warning('Output of applyInverse was size [%s], didn''t match %s sizein: [%s].',...
                    num2str(size(x)),class(this), num2str(this.sizein));
            end
        end
        function M = makeComposition(this, G)
            % Compose the Map \\(\\mathrm{H}\\) with the given Map
            % \\(\\mathrm{G}\\). Returns a new map \\(\\mathrm{M=HG}\\)
            %
            % Calls the method :meth:`makeComposition_`
            if ~cmpSize(this.sizein,G.sizeout)
                error('Input to makeComposition is a %s of sizeout size [%s], which didn''t match the %s sizein [%s].',...
                    class(G),num2str(G.sizeout),class(this), num2str(this.sizein));
            end
            M = this.makeComposition_(G);
        end
        function M = plus(this,G)
            % Overload operator (+) for :class:`Map` objects
            % $$ \\mathrm{M}(\\mathrm{x}) := \\mathrm{H}(\\mathrm{x}) + \\mathrm{G}(\\mathrm{x})$$
            %
            % Calls the method :meth:`plus_`
            if ~cmpSize(this.sizein, G.sizein) % check input size
                error('Input to plus is a %s of sizein  [%s], which didn''t match the added %s sizein [%s].',...
                    class(G),num2str(G.sizein),class(this), num2str(this.sizein));
            end
            if ~cmpSize(this.sizeout, G.sizeout) % check input size
                error('Input to plus is a %s of sizeout  [%s], which didn''t match the added %s sizeout [%s].',...
                    class(G),num2str(G.sizeout),class(this), num2str(this.sizeout));
            end
            M = this.plus_(G);
        end
        function M = minus(this,G)
            % Overload operator (-) for :class:`Map` objects
            % $$ \\mathrm{M}(\\mathrm{x}) := \\mathrm{H}(\\mathrm{x}) - \\mathrm{G}(\\mathrm{x})$$
            %
            % Calls the method :meth:`minus_`
            if ~cmpSize(this.sizein, G.sizein) % check input size
                error('Input to plus is a %s of sizein  [%s], which didn''t match the substracted %s sizein [%s].',...
                    class(G),num2str(G.sizein),class(this), num2str(this.sizein));
            end
            if ~cmpSize(this.sizeout, G.sizeout) % check input size
                error('Input to plus is a %s of sizeout  [%s], which didn''t match the substracted %s sizeout [%s].',...
                    class(G),num2str(G.sizeout),class(this), num2str(this.sizeout));
            end
            M = this.minus_(G);
        end
        function M = mpower(this,p)
            % Returns a new :class:`Map` which is the power p \\(\\mathrm{H}^{p}\\) of the
            % current \\(\\mathrm{H}\\).
            %
            % Calls the method :meth:`mpower_`
            M=this.mpower_(p);
        end
        function G = mtimes(this,G)
            % Overload operator (*) for :class:`Map` objects
            % $$ \\mathrm{M}(\\mathrm{x}) := \\mathrm{H}(\\mathrm{G}(\\mathrm{x}))$$
            %  - If \\(\\mathrm{G}\\) is numeric of size sizein, then :meth:`apply` is called
            %  - If \\(\\mathrm{G}\\) is a :class:`Map`, then a
            %    :class:`MapComposition`is intanciated
            if isa(G,'Map')
                if (isnumeric(this) && isscalar(this)) % Left multiplication by scalar
                    if ~isequal(this,1) % if multiply by 1 do noting
                        if isa(G,'Cost')
                            G=CostMultiplication(this,G);
                        else
                            H=LinOpDiag(G.sizeout,this);
                            G=H.makeComposition(G);
                        end
                    end
                else
                    G =this.makeComposition(G);
                end
            else
                G = this.apply(G);
            end
        end
        function M = times(this,G)
            % Returns a new :class:`Map` which is the element-wise multiplication of the
            % current \\(\\mathrm{H}\\) with \\(\\mathrm{G}\\)
            % $$ \\mathrm{M}(\\mathrm{x}) := \\mathrm{H}(\\mathrm{x}) \\times \\mathrm{G}(\\mathrm{x})$$
            %
            % Calls the method :meth:`times_`
            if ~cmpSize(this.sizein, G.sizein) % check input size
                error('Input to times is a %s of sizein  [%s], which didn''t match the multiplied %s sizein [%s].',...
                    class(G),num2str(G.sizein),class(this), num2str(this.sizein));
            end
            if ~cmpSize(this.sizeout, G.sizeout) % check input size
                error('Input to times is a %s of sizeout  [%s], which didn''t match the multiplied %s sizeout [%s].',...
                    class(G),num2str(G.sizeout),class(this), num2str(this.sizeout));
            end
            M=this.times_(G);
        end
        % TODO Overload operator .^ to do (H(x)).^p ?
        function sz = size(this, varargin)
            sz = {this.sizeout, this.sizein};
            if length(varargin) == 1
                sz = sz{varargin{:}};               
            end
        end
    end
    
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - applyInverse_(this,y)
    % - plus_(this,G)
    % - minus_(this,G)
    % - mpower_(this,p)
    % - times_(this,G)
    % - makeInversion_(this)
    % - makeComposition_(this, H)
    methods (Access = protected)
        function y = apply_(this, x)
            % Not implemented in this Abstract class
            error('apply_ method not implemented');
        end
        function x = applyJacobianT_(this, y, v)
            % Not implemented in this Abstract class
            error('applyJacobianT_ method not implemented');
        end
        function x = applyInverse_(this, y)
            % Not implemented in this Abstract class
            error('applyInverse_ method not implemented');
        end
        function M = plus_(this,G)
            % Constructs a :class:`MapSummation` object to sum the
            % current :class:`Map` \\(\\mathrm{H}\\) with the given \\(\\mathrm{G}\\).
            M = MapSummation({this,G},[1,1]);
        end
        function M = minus_(this,G)
            % Constructs a :class:`MapSummation` object to subtract to the
            % current :class:`Map` \\(\\mathrm{H}\\), the given \\(\\mathrm{G}\\).
            M = this + (-1)*G;
        end
        function M = mpower_(this,p)
            % When \\(p=-1\\), constructs a :class:`MapInversion` object which
            % is the inverse Map of \\(\\mathrm{H}\\).
            % When \\(p\\neq-1\\), this method is not implemented in this Abstract class
            if p==-1
                M=this.makeInversion_();
            else
                error('mpower_ method not implemented for the given power==%d',p);
            end
        end
        function M = times_(this,G)
            % Constructs a :class:`MapMultiplication` object to element-wise multiply the
            % current :class:`Map` \\(\\mathrm{H}\\) with the given \\(\\mathrm{G}\\).
            
            M=MapMultiplication(this,G);
        end
        function M = makeInversion_(this)
            % Constructs a :class:`MapInversion` corresponding to
            % \\(\\mathrm{H}^{-1}\\)
            
            M=MapInversion(this);
        end
        function M = makeComposition_(this, G)
            % Constructs a :class:`MapComposition` object to compose the
            % current Map \\(\\mathrm{H}\\)  with the given \\(\\mathrm{G}\\).
            
            if isa(G,'MapInversion') && isequal(G.M,this)
                M=LinOpScaledIdentity(this.sizein,1);
            elseif isa(G,'MapComposition')
                M=(this*G.H1)*G.H2;
            else
                M = MapComposition(this,G);
            end
        end
    end
    
    %% Utility methods
    % - memoize(this, fieldName, fcn, xs)
    % - clearCaches
    % - clearListenerList
    % - delete
    methods (Access = protected)
        function x = memoize(this, fieldName, fcn, x)
            if ~iscell(x) % handle single input case
                x = {x};
            end
            if   ~isequal(this.memoCache.(fieldName).in, x)
                this.memoCache.(fieldName).in = x;
                x = fcn(x{:});
                this.memoCache.(fieldName).out = x;
            else
                x = this.memoCache.(fieldName).out;
            end
        end
    end
    methods
        function clearCaches(this)
            % Clear precomputation and memoize caches
            this.precomputeCache = struct();       % clean the precompute cache
            fnames = fieldnames(this.memoizeOpts); % clean the memoize cache
            for i = 1:length(fnames)
                this.memoCache.(fnames{i})=struct('in', [], 'out', []);
            end
        end
        function clearListenerList(this)
            % Clear properly the listener list
            prop=properties(this);
            for i=1:length(prop)
                % If public property, add a PostSet listener
                pp=findprop(this,prop{i});
                if ~isempty(pp) && strcmp(pp.SetAccess,'public')
                    if isa(this.(prop{i}),'Map')
                        this.(prop{i}).clearListenerList();
                    elseif isa(this.(prop{i}),'cell') && isa(this.(prop{i}){1},'Map')
                        for j=1:length(this.(prop{i}))
                            this.(prop{i}){j}.clearListenerList();
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
        function delete(this)
            this.clearListenerList();
        end
    end
end