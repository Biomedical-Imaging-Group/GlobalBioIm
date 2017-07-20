classdef LinOpSum <  LinOp
    % Summation linear operator which sums the elements of a variable along
    % given directions.
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{y_k} = \\sum_l \\mathrm{x}_{k,l} $$
    %
    % :param sz: size of \\(\mathrm{x}\\) on which the :class:`LinOpSum` applies.
    % :param index: dimensions along which the sum will be performed (inner sum over l)
    % 
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % See also :class:`LinOp`
    
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
        index  % index along which dimension are computed the finite differences
        ndms   % number of dimensions of the input
        kerdims % ker dimensions
        imdims % im dimensions
    end
    methods
        function this = LinOpSum(sz,index)
            if nargin == 1
                index = [];
            end
            this.name ='LinOp Sum ';
            this.iscomplex= true;
            this.isinvertible=false;
            
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizein = sz;
            
            this.ndms = length(this.sizein);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizein(2) ==1
                this.ndms = 1;
            end
            
            
            if (~isempty(index))
                assert(isvector(index) && length(index)<= this.ndms && max(index)<= this.ndms,'The index should be a conformable  to sz');
                this.index = sort(index,'descend');
            else
                this.index = 1:this.ndms;
            end
            T = true(this.ndms,1);
            T(this.index)=false;
            
            %size of the output = size of the input x length of the index
            % Special case for scalar vectors as matlab thought it is 2D matrix ;-(
            switch(length(this.index))
                case(this.ndms)
                this.sizeout= [1 1];
                case(this.ndms-1)
                this.sizeout= [this.sizein(T) 1];
                otherwise
                this.sizeout= this.sizein(T);
            end
            this.kerdims = this.sizein;
            this.kerdims(T)=1;
            this.imdims = this.sizein;
            this.imdims(~T)=1;
            
        end
        function y = apply(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            
            assert( isequal(size(x),this.sizein),  'x                                                                                                                                                                                   does not have the right size: [%d, %d, %d,%d]',this.sizein);
            for n=this.index
                x = sum(x,n);
            end
            y = squeeze(x);
        end
        function y = adjoint(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            % $$\\mathrm{H}^* : \\mathrm{x} \\mapsto \\mathrm{y_{k,l}} =  \\mathrm{x}_{k} \\; \\forall l$$
            
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
            y = reshape(repmat(reshape(x,this.imdims),this.kerdims),this.sizein);
        end
        function y = HHt(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d]',this.sizeout);
            a = prod(this.kerdims);
            y = x.*a;
        end
    end
end
