classdef LinOp < Map
    % Abstract class for linear operators
	% $$ \\mathrm{H}: \\mathrm{X}\\rightarrow \\mathrm{Y}.$$
    % where \\(\\mathrm{X}\\) and \\(\\mathrm{Y}\\) are either
    % \\(\\mathbb{R}^N\\) or \\(\\mathbb{C}^N\\).
    %
    % All attributes of parent class :class:`Map` are inherited 
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2017 
    %     F. Soulez ferreol.soulez@epfl.ch
    %     M. McCann michael.mccann@epfl.ch  
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
    
   
    %% Constructor
    methods
        function this=LinOp()
            % Add new fields to memoizeOpts and memoCache
            this.memoizeOpts.applyAdjoint=false; 
            this.memoizeOpts.applyAdjointInverse=false; 
            this.memoizeOpts.applyHtH=false;
            this.memoizeOpts.applyHHt=false;
            this.memoCache.applyAdjointInverse=struct('in', [], 'out', []); 
            this.memoCache.applyAdjoint=struct('in', [], 'out', []); 
            this.memoCache.applyHtH=struct('in', [], 'out', []);
            this.memoCache.applyHHt=struct('in', [], 'out', []);
            this.isDifferentiable=true;
        end
    end
   
    
    %% Interface Methods (cannot be overloaded in derived classes: Sealed)
    % In addition to inherited methods from Map
    % - applyAdjoint(this, y)
    % - applyHtH(this,x)
    % - applyHHt(this,y)
    % - transpose(this)
    % - ctranspose(this)
    % - applyAdjointInverse(this,x)
    % - makeHtH(this)
    % - makeHHt(this)
    methods (Sealed)
        function x = applyAdjoint(this, y)
            % Computes  \\(\\mathrm{y=H}^*\\mathrm{y}\\) for \\(\\mathrm{y} \\in \\mathrm{Y}\\)
            %
            % Calls the method :meth:`applyAdjoint_`
            if ~checkSize(y, this.sizeout) % check input size
                error('Input to applyAdjoint was size [%s], didn''t match stated sizeout: [%s].',...
                    num2str(size(y)), num2str(this.sizeout));
            end            
            % memoize
            x = this.memoize('applyAdjoint', @this.applyAdjoint_, y);            
            % check output size
            if ~checkSize(x, this.sizein)
                warning('Output of applyAdjoint was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end
        end
        function y = applyHtH(this,x)
            % Computes  \\(\\mathrm{y=H}^*\\mathrm{Hx}\\) for \\(\\mathrm{y} \\in \\mathrm{Y}\\)
            %
            % Calls the method :meth:`applyHHt_`
            if ~checkSize(x, this.sizein) % check input size
                error('Input to applyHtH was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end            
            % memoize
            y = this.memoize('applyHtH', @this.applyHtH_, x);            
            % check output size
            if ~checkSize(y, this.sizein)
                warning('Output of applyHtH was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(y)), num2str(this.sizein));
            end
        end
        function x = applyHHt(this,y)
            % Computes  \\(\\mathrm{y=HH}^*\\mathrm{y}\\) for \\(\\mathrm{y} \\in \\mathrm{Y}\\)
            %
            % Calls the method :meth:`applyHHt_`
            if ~checkSize(y, this.sizeout) % check input size
                error('Input to applyHHt was size [%s], didn''t match stated sizeout: [%s].',...
                    num2str(size(y)), num2str(this.sizeout));
            end            
            % memoize
            x = this.memoize('applyHHt', @this.applyHHt_, y);            
            % check output size
            if ~checkSize(x, this.sizeout)
                warning('Output of applyHHt was size [%s], didn''t match stated sizeout: [%s].',...
                    num2str(size(x)), num2str(this.sizeout));
            end
        end
        function L = transpose(this)
            % Returns a new :class:`LinOp` which is the Adjoint \\(\\mathrm{H}^{\\star}\\) of the
            % current \\(\\mathrm{H}\\).          
            if this.isComplexOut
                warning('Warning: Do you mean adjoint? For LinOp objects transpose() is an alias of adjoint method');
            end
            L = this.makeAdjoint_();
        end
        function L = ctranspose(this)
            % Do the same as :meth:`transpose`
            L = this.makeAdjoint_();
        end   
        function L = makeAdjoint(this)
            % Do the same as :meth:`transpose`
            L = this.makeAdjoint_();
        end
        function y= applyAdjointInverse(this,x)
            % Computes \\(\\mathrm{y} = \\mathrm{H}^{-\star} \\mathrm{x}\\) for the given
            % \\(\\mathrm{x} \\in \\mathrm{X}\\). (if applicable)
            %
            % Calls the method :meth:`applyAdjointInverse_`
            if ~checkSize(x, this.sizein) % check input size
                error('Input to applyAdjointInverse was size [%s], didn''t match stated sizein: [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end
            % memoize
            y = this.memoize('applyAdjointInverse', @this.applyAdjointInverse_, x);
            % check output size
            if ~checkSize(y, this.sizeout)
                warning('Output of applyAdjointInverse was size [%s], didn''t match stated sizeout: [%s].',...
                    num2str(size(y)), num2str(this.sizeout));
            end
        end
        function M=makeHtH(this)
            % Compose the Adjoint Map \\(\\mathrm{H}^{\\star}\\) with 
            % \\(\\mathrm{H}\\). Returns a new map \\(\\mathrm{M=H}^{\\star} \\mathrm{H}\\)
            %
            % Calls the method :meth:`makeHtH_`
            M=this.makeHtH_();
        end
        function M=makeHHt(this)
            % Compose the  Map \\(\\mathrm{H}\\)  with its adjoint
            % \\(\\mathrm{H}^{\\star}\\). Returns a new map \\(\\mathrm{M=H}\\mathrm{H}^{\\star}\\)
            %
            % Calls the method :meth:`makeHHt_`
            M=this.makeHHt_();
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyAdjoint_(this, y)
    % - applyHtH_(this,x)
    % - applyHHt_(this,y)
    % - applyAdjointInverse_(this,x)
    % - plus_(this,G)
    % - minus_(this,G)
    % - makeAdjoint_(this)
    % - makeHtH_(this)
    % - makeHHt_(this)
    % - makeComposition_(this, G)
    methods (Access = protected) % all the other underscore methods
        function x = applyAdjoint_(this, y)
            % Not implemented in this Abstract class
            error('applyAdjoint_ method not implemented');
        end
        function y = applyHtH_(this,x)
            % There is a default implementation in the abstract class :class:`LinOp` which
            % calls successively the :meth:`apply` and :meth:`applyAdjoint` methods. However, it
            % can be reimplemented in derived classes if there exists a faster way to perform computation.
            y = this.applyAdjoint(this.apply(x));
        end
        function x = applyHHt_(this,y)
            % There is a default implementation in the abstract class :class:`LinOp`
            % which calls successively the :meth:`applyAdjoint` and :meth:`apply` methods.
            % However, it can be reimplemented in derived classes if there exists a faster
            % way to perform computation.
            x = this.apply(this.applyAdjoint(y));
        end
        function y = applyAdjointInverse_(this,x)
            % Not implemented in this Abstract class
            error('applyAdjointInverse_ method not implemented');
        end
        function M = plus_(this,G)
            % If \\(\\mathrm{G}\\) is a :class:`LinOp`, constructs a :class:`LinOpSummation` object to sum the
            % current :class:`LinOp` \\(\\mathrm{H}\\) with the given \\(\\mathrm{G}\\).
            % Otherwise the summation will be a :class:`MapSummation`.
            if isa(G,'LinOp')
                M = LinOpSummation({this,G},[1,1]);
            else
                M = MapSummation({this,G},[1,1]);
            end
        end
        function M = minus_(this,G)
            % If \\(\\mathrm{G}\\) is a :class:`LinOp`, constructs a :class:`LinOpSummation` object to subtract to the
            % current :class:`LinOp` \\(\\mathrm{H}\\), the given \\(\\mathrm{G}\\).
            % Otherwise the summation will be a :class:`MapSummation`.
            if isa(G,'LinOp')
                M = LinOpSummation({this,G},[1,-1]);
            else
                M = MapSummation({this,G},[1,-1]);
            end
        end
        function M = makeAdjoint_(this)
            % Constructs a :class:`LinOpAdjoint` from the current
            % current :class:`LinOp` \\(\\mathrm{H}\\) 
            M=LinOpAdjoint(this);
        end
        function M = makeHtH_(this)
            % Constructs a :class:`LinOpComposition` corresponding to 
            % \\(\\mathrm{H}^{\\star}\\mathrm{H}\\)
            M=LinOpComposition(this',this);
        end
        function M = makeHHt_(this)
            % Constructs a :class:`LinOpComposition` corresponding to 
            % \\(\\mathrm{H}\\mathrm{H}^{\\star}\\)
            M=LinOpComposition(this,this');
        end
        function M = makeComposition_(this, G)
            % Reimplemented from parent class :class:`Map`.
            % If \\(\\mathrm{G}\\) is a :class:`LinOp`, constructs a :class:`LinOpComposition`
            % object to compose the current LinOp (this) with the given :class:`LinOp`\\(\\mathrm{G}\\). 
            % Otherwise the composition will be a :class:`MapComposition`.
            if isa(G,'LinOp')
                if isa(G,'LinOpDiag') && G.isScaledIdentity
                    if G.diag==1
                        M=this;
                    else
                        M=LinOpDiag(this.sizeout,G.diag)*this;  % To always have the diag Op in H1 for LinOpComposition
                    end
                % is HHt or HtH
                elseif (isa(G,'LinOpAdjoint') && isequal(G.TLinOp,this)) || (isa(this,'LinOpAdjoint') && isequal(this.TLinOp,G))
                    M=this.makeHHt();
                % composition with a sum of LinOps
                elseif isa(G,'LinOpSummation')
                    M=G.alpha(1)*this*G.mapsCell{1};
                    for i=2:G.numMaps
                        M=M+G.alpha(i)*this*G.mapsCell{i};
                    end
                % composition with a LinOpComposition
                elseif isa(G,'LinOpComposition') 
                    % to handle properly scalar multiplications
                    if isa(G.H1,'LinOpDiag') && G.H1.isScaledIdentity
                        if isa(this,'LinOpDiag') && this.isScaledIdentity
                            M=LinOpDiag(this.sizeout,G.H1.diag*this.diag)*G.H2;
                        else
                            M=LinOpDiag(this.sizeout,G.H1.diag)*(this*G.H2);
                        end
                    else
                        M=this*G.H1*G.H2;
                    end
                % composition with other LinOps
                else
                    M = LinOpComposition(this,G);
                end
            else
                M = makeComposition_@Map(this,G);
            end
        end
        function M = mpower_(this,p)
            % Reimplemented from :class:`Map`  
            if p==-1
                M=LinOpInversion(this);
            else
                M=mpower_@Map(this,p);
            end
        end
    end
    
    %% Methods of superclass Map that do not need to be reimplemented in derived Costs
    % - applyJacobianT_(this, y, v)
    methods (Access = protected, Sealed)
        function x = applyJacobianT_(this, y, v)
            % Uses the method applyAdjoint (hence do not need to be
            % reimplemented in derived classes)
            x = this.applyAdjoint(y);
        end
    end
          
		% TODO: ask ferreol: is this needed? when do we have fast HT W H?
%         function y = HtWH(this,x,W) 
%         	% Apply \\(\\mathrm{H}^*\\mathrm{WH}\\)
%         	%
%         	% :param x: \\(\\in X\\)
%         	% :param W: a :class:`LinOp` object
%         	% :returns y: \\(= \\mathrm{H^*WHx}\\)
%         	%
% 
%             if (isscalar(W) && isreal(W))
%                 y = W.*this.HtH(x);
%             else
%                 assert(isa(W,'LinOp'),'W must be a LinOp');
%                 y = this.adjoint(W.apply(this.apply(x)));
%             end
%         end           
end

