classdef OneToMany < LinOp
    %% OneToMany :  linear operators
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = OneToMany(ALinOp,alpha)
    % Element wise sum  of LinOps:
    % Array of LinOp that will be applied to the same vector x contained in cell array ALINOP weighted by ALPHA
    % (default 1)
    % Obj{n} = alpha(n) * ALinOp{n}
    % Applying this operator returning the cell array:
    % Obj.Apply(x) = alpha{n} * Obj.ALinOp{n}.Apply(x)
    % WARNING
    % Obj.HtH(x) = alpha{n} * Obj.ALinOp{n}.HtH(x)
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
        ALinOp     % Cell array of linop
        numLinOp   % number of linop
        alpha      % scalar factor
    end
    
    methods
        function this = OneToMany(ALinOp,alpha)
            this.name ='OneToMany';
            if nargin == 1
                alpha = 1;
            end
            
            
            assert(iscell(ALinOp) && isa(ALinOp{1}(1),'LinOp'),'First input should be a cell array LinOp');
            this.ALinOp = ALinOp;
            this.numLinOp = numel(ALinOp);
            assert(isnumeric(alpha)&& ( isscalar(alpha) || ( isvector(alpha) && (numel(alpha)== this.numLinOp))),'second input should be a scalar or an array of scalar of the same size as the first input');
            if  isscalar(alpha)
                this.alpha = repmat(alpha,this.numLinOp,1) ;
            else
                this.alpha = alpha;
            end
            this.iscomplex= this.ALinOp{1}(1).iscomplex;
            this.issquare = this.ALinOp{1}(1).issquare;
            this.isinvertible=false;
            this.sizein =  this.ALinOp{1}(1).sizein;
            for n =2:this.numLinOp
                assert(isempty(this.ALinOp{n}(1).sizein) ||isequal(this.sizein,this.ALinOp{n}(1).sizein),'%d-th input does not have the right right hand side size ') ;
                this.iscomplex= this.ALinOp{n}(1).iscomplex ||  this.iscomplex ;
            end
            
            
            
        end
        
        function y = Apply(this,x) % Apply the operator
               assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);
            y = cell(1,this.numLinOp);
            for n =1:this.numLinOp
                y{n} = this.alpha(n) .* this.ALinOp{n}(1).Apply(x);
            end
        end
        function y = Adjoint(this,x) % Apply the adjoint
            assert( iscell(x), 'Adjoint of a OneToMany must be applied to a cell array',this.sizeout);
           assert( numel(x)==this.numLinOp, 'OneToMany LinOp: input does not have the right number of cell: %d',this.numLinOp); 
            y =  this.alpha(1) .* this.ALinOp{1}(1).Adjoint(x{1});
            for n =2:this.numLinOp
            assert( isequal(size(x{n}),this.ALinOp{n}.sizeout),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizeout);
                y = y + this.alpha(n) .* this.ALinOp{n}(1).Adjoint(x{n});
            end
        end
        function y = HtH(this,x) % Apply the operator
               assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);
            y = this.alpha(1) .* this.ALinOp{1}(1).HtH(x);
            for n =2:this.numLinOp
                y = y + this.alpha(n) .* this.ALinOp{n}(1).HtH(x);
            end
        end
        function y = HtWH(this,x,W) % Apply the operator 
           
            assert( iscell(W), 'W in OneToMany.HtWH must be a cell array',this.sizeout);
           assert( numel(W)==this.numLinOp, 'OneToMany W does not have the right number of cell: %d',this.numLinOp);
               assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d %d %d %d ]',this.sizein);
            y = this.alpha(1) .* this.ALinOp{1}(1).HtWH(x,W{1});
            for n =2:this.numLinOp
                y = y + this.alpha(n) .* this.ALinOp{n}(1).HtWH(x,W{n});
            end
        end
    end
end

