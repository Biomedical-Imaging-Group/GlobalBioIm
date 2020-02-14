classdef LinOpBroadcastMatrix <  LinOp
    % LinOpBroadCastMatrix which applies a matrix along a given dimension
    % 
    % Build a LinOp from a given matrix 
    %
    % :param M: of size M x N
    % :param sz: input size of the :class:`LinOpBroadCastMatrix`
    % :param index: dimension along which the matrix is applied (sz(index) must be equal to N)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % If M is 2x3, sz=[2 3 2] and index=2, then the matrix M is applied to
    % each x(i,:,j) leading to an output of size [2 2 2].
    %
    % **Example** H=LinOpBroadcastMatrix(M,sz,index)
    %
    % See also :class:`LinOp`, :class:`Map`
     
    %%    Copyright (C) 2018 
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
    
    properties (SetAccess = protected,GetAccess = public)
        M            % matrix
        index        % index along which the matrix is applied
        indexdiff;   % other idexes
        ndms         % number of dimensions
        permIdxIn;   % indexes for input permutation
        permIdxOut;  % indexes for input permutation
    end
    
    %% Constructor
    methods
        function this = LinOpBroadcastMatrix(M,sz,index)
            this.name='LinOpBroadcastMatrix';
            this.M=M;  
            if nargin<2
                sz=[size(M,2),1];
                index=1;
            elseif nargin <3
                index=find(sz~=1,1,'last');
            end
            this.ndms=length(sz);
            assert(isscalar(index) && index <= length(sz),'index must be a scalar smaller than the length of sz');
            assert(size(M,2)==sz(index),'size(M,2) and sz(index) must be equal');
            this.index=index;
            this.indexdiff=setdiff(1:this.ndms,index);
            this.sizein=sz;
            this.sizeout=sz;this.sizeout(index)=size(M,1);
            this.permIdxIn=[index,this.indexdiff];
            this.permIdxOut=[2:index,1,index+1:this.ndms];
            this.norm = norm(M);                      
            this.isDifferentiable=true;
            if size(M,1)==size(M,2) && det(M)~=0
                this.isInvertible=true;
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function x=apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            
            x=reshape(permute(x,this.permIdxIn),[this.sizein(this.index),prod(this.sizein(this.indexdiff))]);
            x=permute(reshape(this.M*x,[this.sizeout(this.index),this.sizeout(this.indexdiff)]),this.permIdxOut);
        end       
        function y=applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            
            tmp=reshape(permute(x,this.permIdxIn),[this.sizeout(this.index),prod(this.sizeout(this.indexdiff))]);
            y=permute(reshape(this.M'*tmp,[this.sizein(this.index),this.sizein(this.indexdiff)]),this.permIdxOut);
        end
        function y=applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            
            if this.isInvertible
                tmp=reshape(permute(x,this.permIdxIn),[this.sizeout(this.index),prod(this.sizeout(this.indexdiff))]);
                y=permute(reshape(this.M\tmp,[this.sizein(this.index),this.sizein(this.indexdiff)]),this.permIdxOut);
            else
                y=applyInverse_@LinOp(this,x);
            end
        end
        function y=applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            
             tmp=reshape(permute(x,this.permIdxIn),[this.sizein(this.index),prod(this.sizein(this.indexdiff))]);
             y=permute(reshape(this.M'*this.M*tmp,[this.sizein(this.index),this.sizein(this.indexdiff)]),this.permIdxOut);     
        end
        function y=applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            
            tmp=reshape(permute(x,this.permIdxIn),[this.sizeout(this.index),prod(this.sizeout(this.indexdiff))]);
            y=permute(reshape(this.M*this.M'*tmp,[this.sizeout(this.index),this.sizeout(this.indexdiff)]),this.permIdxOut);
        end
        function H=makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOp`.
            
            H=LinOpBroadcastMatrix(this.M',this.sizeout,this.index);
        end
        function H=makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            
            H=LinOpBroadcastMatrix(this.M'*this.M,this.sizein,this.index);
        end
        function H=makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            
            H=LinOpBroadcastMatrix(this.M*this.M',this.sizeout,this.index);
        end
        function H = makeInversion_(this)
            % Reimplemented from parent class :class:`LinOp`.
            
            if this.isInvertible
                H=LinOpBroadcastMatrix(inv(this.M),this.sizeout,this.index);
            else
                H=makeInversion_@LinOp(this);
            end
        end
        function H = plus_(this,G)
            % Reimplemented from parent class :class:`LinOp`.
            
            if isa(G,'LinOpBroadcastMatrix') && this.index==G.index
                H=LinOpBroadcastMatrix(this.M+G.M,this.sizein,this.index);
            elseif isa(G,'LinOpDiag') && G.isScaledIdentity
                H=LinOpBroadcastMatrix(this.M+G.diag*eye(size(this.M)),this.sizein,this.index);
            else
                H=plus_@LinOp(this,G);
            end
        end
        function H = makeComposition_(this, G)
            % Reimplemented from parent class :class:`LinOp`.
                   
            if isa(G,'LinOpBroadcastMatrix') && this.index==G.index && size(this.M,2)==size(G.M,1)
                H=LinOpBroadcastMatrix(this.M*G.M,G.sizein,G.index);
            elseif isa(G,'LinOpDiag') && G.isScaledIdentity
                H=LinOpBroadcastMatrix(this.M*G.diag,this.sizein,this.index);
            else
                H=makeComposition_@LinOp(this,G);
            end
        end
    end
end

