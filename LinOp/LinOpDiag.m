classdef LinOpDiag <  LinOp
    % LinOpDiag: Diagonal operator
    % $$ \\mathrm{Hx}= \\mathrm{\\mathbf{diag}(w)x}$$
    % where \\(\\mathrm{w} \\in \\mathbb{R}^N\\) or \\(\\mathbb{C}^N\\) is a
    % vector containing the diagonal elements of \\(\\mathrm{H}\\).
    %
    % :param diag: elements of the diagonal (Non-singleton dimensions of diag and sz must be consistent)
    % :param sz: size (if not the size of the given diag).
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
    
    %% Properties
    % - Public
    properties (SetObservable, AbortSet)
        diag;                    % diagonal or a scalar
    end
    % - Readable
    properties (SetAccess = protected,GetAccess = public)
        isScaledIdentity=false;  % if diag is constant, it is not stored
    end
    
    %% Constructor
    methods
        function this = LinOpDiag(sz,diag)
            % Default values
            if nargin <1, error('At least a size should be given');end
            if nargin <2, diag=1; end
            if isempty(sz), sz=size(diag); end
            % Checks
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            % Set properties
            this.name ='LinOpDiag';
            this.sizeout=sz;
            this.sizein=sz;         
            this.diag=diag;
            this.isDifferentiable=true;
            % Initialize 
            this.initObject('LinOpDiag');
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`LinOp`
            
            % Call superclass method
            updateProp@LinOp(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'diag') || strcmp(prop,'all')
                assert(isnumeric(this.diag),'diag must be numeric');
                assert(length(this.sizein) >= length(size(this.diag)),'Number of dimensions of diag must be smaller than the one of the given size');
                assert(all(1-((this.sizein(1:length(size(this.diag)))-size(this.diag)).*(1-(size(this.diag)==1)))),'Non-singleton dimensions of diag and sizein must match each other.');     
                if sum(this.diag(:)==0)==0
                    this.isInvertible=true;
                else
                    this.isInvertible=false;
                end
                if isscalar(this.diag) || norm(this.diag(:)-this.diag(1))<1e-13
                    this.isScaledIdentity=true;
                    this.diag = this.diag(1);
                end
                this.norm=max(abs(this.diag(:)));
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function x = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.            
                if  verLessThan('matlab', '9.1')
                    x =bsxfun(@times,this.diag,x);
                else
                    x=this.diag.*x;
                end
        end
        function x = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
             if  verLessThan('matlab', '9.1')
                 x =bsxfun(@times,conj(this.diag),x);
             else
                 x=conj(this.diag).*x;
             end
        end
        function x = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if verLessThan('matlab', '9.1')
                x =bsxfun(@times,abs(this.diag).^2,x);
            else
                x=abs(this.diag).^2.*x;                
            end
        end
        function x = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            x=this.applyHtH(x);
        end
        function x = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                if verLessThan('matlab', '9.1') 
                    x =bsxfun(@times,(1./this.diag),x);                    
                else
                    x =x./this.diag;                   
                end
            else
                x = applyInverse_@LinOp(this,x);
            end
        end
        function x = applyAdjointInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                if verLessThan('matlab', '9.1')
                    x =bsxfun(@times,(1./conj(this.diag)),x);
                else
                    x =x./conj(this.diag);
                end
            else
                x = applyAdjointInverse_@LinOp(this,x);
            end
        end
        function M = plus_(this,G)
            % Reimplemented from parent class :class:`LinOp`.
            if isa(G,'LinOpDiag')
                M=LinOpDiag(this.sizein,bsxfun(@plus,G.diag,this.diag));
            elseif isa(G,'LinOpConv') && this.isScaledIdentity
                if G.useRFT
                    M=LinOpConv(this.diag+G.mtf,G.isReal,G.index,'useRFT');
                else
                    M=LinOpConv(this.diag+G.mtf,G.isReal,G.index);
                end
            elseif isa(G,'LinOpMatrix') && this.isScaledIdentity
                M=LinOpMatrix(G.M+this.diag*eye(size(G.M)),G.sizein,G.index);
            else
                M=plus_@LinOp(this,G);
            end
        end
% Not needed anymore : will goes at the level of Map and use sum with a wight (-1).     
%         function M = minus_(this,G)
%             % Reimplemented from parent class :class:`LinOp`.
%             if isa(G,'LinOpDiag')
%                 M=LinOpDiag(this.sizein,bsxfun(@minus,this.diag,G.diag));
%             elseif isa(G,'LinOpConv') && this.isScaledIdentity
%                 if G.useRFT
%                     M=LinOpConv(this.diag-G.mtf,G.isReal,G.index,'useRFT');
%                 else
%                     M=LinOpConv(this.diag-G.mtf,G.isReal,G.index);
%                 end
%             else
%                 M=minus_@LinOp(this,G);
%             end
%         end
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
            
            sz=ones(size(this.sizein));   % to consider singleton dimension that are at the end
            sz(1:length(size(this.diag)))=size(this.diag);
            if isa(G,'LinOpDiag')
                M=LinOpDiag(this.sizein,G.diag.*this.diag);
            elseif isa(G,'LinOpConv') && ( this.isScaledIdentity || (all(sz(G.index)==1)) )
                if G.useRFT
                    M = LinOpConv('PSF',iSrft(this*G.mtf,G.Notindex),G.isReal,G.index,'useRFT');
                else
                    M = LinOpConv(this*G.mtf,G.isReal,G.index);
                end           
            elseif isa(G,'LinOpBroadcastMatrix') && this.isScaledIdentity
                M=LinOpBroadcastMatrix(G.M*this.diag,G.sizein,G.index);
            else
                M=makeComposition_@LinOp(this,G);
            end
        end
    end
end
