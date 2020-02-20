classdef StackLinOp < LinOp
    %% StackLinOp : Stack of linear operators
    %  Matlab Linear Operator Library
    %
    % TODO: INTEGRATE THIS IN THE NEW VERSION
    %
    % Example
    % Obj = StackLinOp(ALinOp,alpha)
    % Stack all linop contained in vector ALINOP weighted by ALPHA
    % (default 1)
    % such that
    % y(..,i) = alpha(i) * ALinOp{i}
    %
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
        ALinOp     % Array of linop
        numLinOp   % number of linop
        alpha      % scalar factor
        usecomplex = true; % if false complex are represented as an extra dimension of size 2 containning Real and imagenary parts of x
        prmtIndex;
    end
    
    methods
        function this = StackLinOp(ALinOp,alpha, varargin)
            this.name ='StackLinOp';
            
            if nargin == 1
                alpha = 1;
            end
            
            this.numLinOp = numel(ALinOp);
            assert(isnumeric(alpha)&& ( isscalar(alpha) || ( isvector(alpha) && (numel(alpha)== this.numLinOp))),'second input should be a scalar or an array of scalar of the same size as the first input');
            if  isscalar(alpha)
                this.alpha = repmat(alpha, 1, this.numLinOp) ;
            else
                this.alpha = alpha;
            end
            
            allLinOps = all( cellfun(@(x)(isa(x, 'LinOp')), ALinOp) );
            assert(iscell(ALinOp) && allLinOps, 'First input should be a cell array LinOp');
            
            
            this.ALinOp = ALinOp;
            this.isComplex= this.ALinOp{1}(1).isComplex;
            this.isInvertible=false;
            this.sizein =  this.ALinOp{1}(1).sizein;
            this.sizeout =  [this.ALinOp{1}(1).sizeout this.numLinOp];
            for n =2:this.numLinOp
                assert(isempty(this.ALinOp{n}(1).sizein)  || cmpSize(this.sizein,this.ALinOp{n}(1).sizein),'%d-th input does not have the right hand side size ', n) ;
                assert(isempty(this.ALinOp{n}(1).sizeout) ||cmpSize(this.ALinOp{1}(1).sizeout,this.ALinOp{n}(1).sizeout),'%d-th input does not have the left hand side size ', n);
                this.isComplex= this.ALinOp{n}(1).isComplex || this.isComplex ;
            end
            
            for c=1:length(varargin)
                switch varargin{c}
                    case('DontUseComplex')
                        assert(this.sizein(end)==2, ' last dimension of input LinOp should be 2');
                        this.usecomplex = false;
                        this.sizeout = [this.sizeout(1:end-2),this.sizeout(end),this.sizeout(end-1)];
                        this.isComplex= false;
                        nd= numel(this.sizeout);
                        this.prmtIndex = [1:(nd-2) nd (nd-1)];
                end
            end
            
            
            
        end
    end
    methods (Access = protected)
        function y = apply_(this,x) % apply the operator
            assert(checkSize(x, this.sizein));
            
            y = zeros_(prod(this.ALinOp{1}.sizeout),this.numLinOp);
            for n = 1:this.numLinOp
                tmp =this.ALinOp{n}(1).apply(x);
                y(:,n) =  this.alpha(n) .* tmp(:);
            end
            
            if ~this.usecomplex
                y = reshape(y, [this.ALinOp{1}.sizeout this.numLinOp]);
                y = permute(y,this.prmtIndex);
            else
                y = reshape(y, this.sizeout);
            end
        end
        function y = applyAdjoint_(this,x) % apply the adjoint
            assert(checkSize(x, this.sizeout));
            y =  zeros_(this.sizein);
            
            if ~this.usecomplex
                x = permute(x,this.prmtIndex);
            end
            x = reshape(x, prod(this.ALinOp{1}.sizeout),this.numLinOp);
            for n = 1:this.numLinOp
                xtmp = zeros_(this.ALinOp{1}.sizeout);
                xtmp(:) = x(:,n);
                y = y + this.alpha(n) .* this.ALinOp{n}(1).applyAdjoint(xtmp);
            end
        end
    end
    
    methods (Access = protected)
        %% Copy
      function this = copyElement(obj)
          this = copyElement@LinOp(obj);
            for n = 1:this.numLinOp
                 this.ALinOp{n} = copy(obj.ALinOp{n});
            end
      end
    end
end

