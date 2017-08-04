classdef (Abstract) Map < handle
    % Abstract class for Maps which maps elements from \\(\\mathrm{X}\\) to
    % \\(\\mathrm{Y}\\)
    % $$ \\mathrm{H}: \\mathrm{X}\\rightarrow \\mathrm{Y}.$$
    %
    % :param name: name of the linear operator \\(\\mathbf{H}\\)
    % :param sizein:  dimension of the left hand side vector space \\(\\mathrm{X}\\)
    % :param sizeout:  dimension of the right hand side vector space \\(\\mathrm{Y}\\)
    % :param norm: norm of the operator \\(\\|\\mathrm{H}\\|\\) (if known, otherwise -1)
    % :param memoizeOpts: 
    % :param doPrecomputation: boolean true to allow doing precomputations to save 
    % time (will generally require more memory).
     
    
    %     Copyright (C) 2017 M. McCann michael.mccann@epfl.ch & 
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
        name = 'none'           % name of the linear operator
        
        isInvertible = false;     % true if H.applyInverse(  ) will work %todo fix capitalization everywhere
        isDifferentiable = false; % true if H.applyJacobianT(   ) will work
        
        % todo: decide on how to handle checking for inverse, grad, ... being
        % implemented and working (no divide by zero).
        %implementedMetods = struct('applyJacobianT', false);
        
        isComplex = false;      % true is the operator is complex %TODO fix capitalization everywhere
        
        norm=-1;                % norm of the operator
        
        sizein;                 % dimension of the right hand side vector space
        sizeout;                % dimension of the left hand side vector space
        
        memoizeOpts = struct('apply', false);
        doPrecomputation = false;
    end
    
    properties (SetAccess = private)
        memoCache = struct('apply', struct('in', [], 'out', []));
        precomputeCache = struct();
    end
    
    %% Interface Methods (cannot be overloaded in derived classes: Sealed)
    % - apply(this,x)
    % - applyJacobianT(this, y, v)
    % - applyInverse(this,y)
    % - makeComposition(this, G)
    methods (Sealed)
        function y = apply(this, x)
            % Computes \\(\\mathrm{y}=\\mathrm{H}(\\mathrm{x})\\) for the given
            % \\(\\mathrm{x} \\in \\mathrm{X}\\).
            
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
            % \\(\\mathrm{y} \\in \\mathrm{X}\\). (if applicable)
            
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
            % \\(\\mathrm{G}\\). Return a new map \\(\\mathrm{M=HG}\\)
            assert(isa(G,'Map'),'Composition is only with a Map');
            assert(isequal(G.sizeout,this.sizein),'The Map G to compose has sizeout inconsistent with the sizein of this.H');
            M = this.makeComposition_(this, G);          
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - applyInverse_(this,y)
    % - makeComposition_(this, H)
    methods (Access = protected)
        function y = apply_(this, x)
            error('apply_ method not implemented');
        end
        
        function x = applyJacobianT_(this, y, v)
            error('applyJacobianT_ method not implemented');
        end
        
        function x = applyInverse_(this, y)
            error('applyInverse_ method not implemented');
        end
        
        function M = makeComposition_(this, G)
            M = MapComposition({this, G});
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