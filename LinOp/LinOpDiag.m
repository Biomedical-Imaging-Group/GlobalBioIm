classdef LinOpDiag <  LinOp
    % LinOpDiag: Diagonal operator
    % $$ \\mathrm{Hx}= \\mathrm{\\mathbf{diag}(w)x}$$
    % where \\(\\mathrm{w} \\in \\mathbb{R}^N\\) or \\(\\mathbb{C}^N\\) is a
    % vector containing the diagonal elements of \\(\\mathrm{H}\\).
    %
    % :param diag: elements of the diagonal (Non-singleton dimensions of diag and sz must be consistent)
    % :param sz: size (if not the size of the given diag) to build a scaled identity operator.
    %
    % All attributes of parent class :class:`LinOp` are inherited.
    %
    % If size(diag)=[2 1] and sz=[2 2], then LinOpDiag multiplies each
    % column of a given 2x2 matrix with the diag vector (see bsxfun).
    %
    % **Example** D=LinOpDiag(sz,diag)
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %%    Copyright (C) 2015
    %     F. Soulez ferreol.soulez@epfl.ch
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
        diag;                    % diagonal or a scalar
        isScaledIdentity=false;  % if diag is constant, it is not stored
    end
    
    %% Constructor
    methods
        function this = LinOpDiag(sz,diag)
            this.name ='LinOpDiag';
            if nargin <1, error('At least a size should be given');end
            if nargin <2, diag=1; end
            if isempty(sz), sz=size(diag); end
            assert(isnumeric(diag),'diag must be numeric');
            assert(length(sz) >= length(size(diag)),'Number of dimensions of diag must be smaller than the one of the given size');
            assert(all(1-((sz(1:length(size(diag)))-size(diag)).*(1-(size(diag)==1)))),'Non-singleton dimensions of diag and sz must match each other.');
            this.sizeout=sz;
            this.sizein=sz;
            this.isDifferentiable=true;
            if all(diag)
                this.isInvertible=true;
            else
                this.isInvertible=false;
            end
            if isscalar(diag) || norm(diag(:)-diag(1))<1e-13
                this.isScaledIdentity=true;
                this.diag = diag(1);
            else
                this.diag=diag;
            end
            % -- Norm of the operator
            this.norm=max(abs(diag(:)));
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =bsxfun(@times,this.diag,x);
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =bsxfun(@times,conj(this.diag),x);
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =bsxfun(@times,abs(this.diag).^2,x);
        end
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=this.applyHtH(x);
        end
        function y = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                y =bsxfun(@times,(1./this.diag),x);
            else
                y = applyInverse_@LinOp(this,x);
            end
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                y =bsxfun(@times,conj(1./this.diag),x);
            else
                y = applyAdjointInverse_@LinOp(this,x);
            end
        end
        function M = plus_(this,G)
            % Reimplemented from parent class :class:`LinOp`.
            if isa(G,'LinOpDiag')
                M=LinOpDiag(this.sizein,bsxfun(@plus,G.diag,this.diag));
            elseif isa(G,'LinOpConv') && this.isScaledIdentity
                M=LinOpConv(this.diag+G.mtf,G.isReal,G.index);
            else
                M=plus_@LinOp(this,G);
            end
        end
        function M = minus_(this,G)
            % Reimplemented from parent class :class:`LinOp`.
            if isa(G,'LinOpDiag')
                M=LinOpDiag(this.sizein,bsxfun(@minus,this.diag,G.diag));
            elseif isa(G,'LinOpConv') && this.isScaledIdentity
                M=LinOpConv(this.diag-G.mtf,G.isReal,G.index);
            else
                M=minus_@LinOp(this,G);
            end
        end
        function M = makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein,conj(this.diag));
        end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein,abs(this.diag).^2);
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein,abs(this.diag).^2);
        end
        function M = mpower_(this,p)
            % Reimplemented from parent class :class:`LinOp`
            if p>=0 || this.isInvertible
                M=LinOpDiag(this.sizein,(this.diag).^p);
            else
                M=mpower_@LinOp(this,p);
            end
        end
        function  M = makeInversion_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                M=LinOpDiag(this.sizein,1./this.diag);
            end
        end
        function M = makeComposition_(this, G)
            % Reimplemented from parent class :class:`LinOp`.
            
            sz=size(this.diag);
            if isa(G,'LinOpDiag')
                M=LinOpDiag(this.sizein,G.diag.*this.diag);
            elseif isa(G,'LinOpConv') && ( this.isScaledIdentity || (~isempty(G.index) && all(sz(G.index)==1)) )
                M = LinOpConv(this*G.mtf,G.isReal,G.index);
            else
                M=makeComposition_@LinOp(this,G);
            end
        end
    end
end
