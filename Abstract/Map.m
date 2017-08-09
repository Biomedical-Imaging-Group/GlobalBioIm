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
    % :param isComplexIn:  true if \\(\\mathrm{X}\\) is complex valued
    % :param isComplexOut:  true if \\(\\mathrm{Y}\\) is complex valued
    % :param memoizeOpts: 
    % :param doPrecomputation: boolean true to allow doing precomputations to save 
    % time (will generally require more memory).
    
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
    
    %% Public properties
    properties         
        name = 'none'             % name of the linear operator        
        isInvertible = false;     % true if H.applyInverse(  ) will work 
        isDifferentiable = false; % true if H.applyJacobianT(   ) will work
        isComplexIn = false;      % true is the space X is complex valued
        isComplexOut = false;     % true is the space Y is complex valued
        sizein;                   % dimension of the right hand side vector space
        sizeout;                  % dimension of the left hand side vector space   
        norm=-1;                  % norm of the operator    
        memoizeOpts = struct('apply', false, ...
                             'applyJacobianT', false, ...
                             'applyInverse', false);
        doPrecomputation = false;
    end
    
    %% Private properties
    properties (SetAccess = protected,GetAccess = protected)
        memoCache = struct('apply', struct('in', [], 'out', []),...
                           'applyJacobianT', struct('in', [], 'out', []), ...
                           'applyInverse', struct('in', [], 'out', []));
        precomputeCache = struct();
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
	% - size(this, [dim])
    methods (Sealed)
        function y = apply(this, x)
            % Computes \\(\\mathrm{y}=\\mathrm{H}(\\mathrm{x})\\) for the given
            % \\(\\mathrm{x} \\in \\mathrm{X}\\). 
            %
            % Calls the method :meth:`apply_`
            
            if ~checkSize(x, this.sizein) % check input size
                error('Input to apply was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end            
            % memoize
            y = this.memoize('apply', @this.apply_, x);            
            % check output size
            if ~checkSize(y, this.sizeout)
                warning('Output of apply was size [%s], didn''t match stated sizeout: [%s].',...
                    num2str(size(y)), num2str(this.sizeout));
            end
        end   
        function x = applyJacobianT(this, y, v)
            % Compute \\(\\mathrm{x}=[\\mathrm{J}_{\\mathrm{H}}(\\mathrm{v})]^{\\star}\\mathrm{y}\\)
            % where
            %
            % - \\([\\mathrm{J}_{\\mathrm{H}}(\\mathrm{v})]\\) is the Jacobian matrix of
            %   the Map \\(\\mathrm{H}\\) computed at \\(\\mathrm{v} \\in \\mathrm{X} \\)
            % - \\(\\mathrm{y} \\in \\mathrm{Y} \\)
            %
            % Calls the method :meth:`applyJacobianT_`
            
            if ~checkSize(y, this.sizeout) % check input size
                error('Input of y applyJacobianT was size [%s], didn''t match stated sizeout: [%s].',...
                    num2str(size(y)), num2str(this.sizeout));
            end
            if ~checkSize(v, this.sizein)
                error('Input to v applyJacobianT was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(v)), num2str(this.sizein));
            end            
            % memoize
            x = this.memoize('applyJacobianT', @this.applyJacobianT_, {y, v});            
            % check output size
            if ~checkSize(x, this.sizein)
                warning('Output of applyJacobianT was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end
        end    
        function x = applyInverse(this, y)
            % Computes \\(\\mathrm{x} = \\mathrm{H}^{-1} \\mathrm{y}\\) for the given
            % \\(\\mathrm{y} \\in \\mathrm{Y}\\). (if applicable)   
            %
            % Calls the method :meth:`applyInverse_`
            if ~checkSize(y, this.sizeout) % check input size
                error('Input to applyInverse was size [%s], didn''t match stated sizeout: [%s].',...
                    num2str(size(y)), num2str(this.sizeout));
            end           
            % memoize
            x = this.memoize('applyInverse', @this.applyInverse_,y);           
            % check output size
            if ~checkSize(x, this.sizein)
                warning('Output of applyInverse was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end
        end      
        function M = makeComposition(this, G)
            % Compose the Map \\(\\mathrm{H}\\) with the given Map
            % \\(\\mathrm{G}\\). Returns a new map \\(\\mathrm{M=HG}\\)
            %
            % Calls the method :meth:`makeComposition_`
            M = this.makeComposition_(G);          
        end
        function M = plus(this,G)
            % Overload operator (+) for :class:`Map` objects
            % $$ \\mathrm{M}(\\mathrm{x}) := \\mathrm{H}(\\mathrm{x}) + \\mathrm{G}(\\mathrm{x})$$
            %
            % Calls the method :meth:`plus_`
            M = this.plus_(G);
        end       
        function M = minus(this,G)
            % Overload operator (-) for :class:`Map` objects
            % $$ \\mathrm{M}(\\mathrm{x}) := \\mathrm{H}(\\mathrm{x}) - \\mathrm{G}(\\mathrm{x})$$  
            %
            % Calls the method :meth:`minus_`
            M = this.minus_(G);
        end      
        function M = mpower(this,p)
        end
        function M = mtimes(this,G)
            % Overload operator (*) for :class:`Map` objects
            % $$ \\mathrm{M}(\\mathrm{x}) := \\mathrm{H}(\\mathrm{G}(\\mathrm{x}))$$
            %  - If \\(\\mathrm{G}\\) is numeric of size sizein, then :meth:`apply`is called
            %  - If \\(\\mathrm{G}\\) is a :class:`Map`, then a
            %    :class:`MapComposition`is intanciated
            if isa(G,'Map')
                if (isnumeric(this) && isscalar(this)) % Left multiplication by scalar
                    H=LinOpScaledIdentity(G.sizeout,this);
                    M=H.makeComposition(G);
                else
                    M =this.makeComposition(G);
                end
            else
                M = this.apply(G);
            end
        end
        % TODO Overload operator .* to multiply term by term two Maps ? (avoid the use of LinOpDiag)
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
    % - makeComposition_(this, H)
    % - plus_(this,G)
    % - minus_(this,G)
    % - mpower_(this,p)
    % - mtimes(this,G)
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
        function M = makeComposition_(this, G)
            % Constructs a :class:`MapComposition` object to compose the
            % current Map \\(\\mathrm{H}\\)  with the given \\(\\mathrm{G}\\). 
            M = MapComposition(this,G);
        end
        function M = plus_(this,G)
            % Constructs a :class:`MapSummation` object to sum the
            % current Map \\(\\mathrm{H}\\) with the given \\(\\mathrm{G}\\). 
            M = MapSummation({this,G},[1,1]);
        end
        function M = minus_(this,G)
            % Constructs a :class:`MapSummation` object to subtract to the
            % current Map \\(\\mathrm{H}\\), the given \\(\\mathrm{G}\\). 
            M = MapSummation({this,G},[1,-1]);
        end
        function M = mpower_(this,p)
            % When \\(p=-1\\), constructs a :class:`MapInversion` object which
            % is the inverse Map of \\(\\mathrm{H}\\).
            % When \\(p\\neq-1\\), this method is not implemented in this Abstract class
            if p==-1
                M=MapInversion(this);
            else
                error('mpower_ method not implemented');
            end
		end
	end

		
    %% Utility methods  
    % - memoize(this, fieldName, fcn, xs)
    methods (Access = protected)          
        function y = memoize(this, fieldName, fcn, xs)           
            if ~iscell(xs) % handle single input case
                xs = {xs};
            end          
            if ~this.memoizeOpts.(fieldName) || ~isequal(this.memoCache.(fieldName).in, xs)
                y = fcn(xs{:});
                if this.memoizeOpts.(fieldName)
                    this.memoCache.(fieldName).in = xs;
                    this.memoCache.(fieldName).out = y;
                end
            else
                y = this.memoCache.(fieldName).out;
            end
        end    
    end 
end