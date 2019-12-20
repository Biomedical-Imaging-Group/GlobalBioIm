classdef LinOpBroadcast <  LinOp
    % LinOpSum: Broacasting linear operator which broadcast the elements of a variable along
    % given directions.
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{y_{k,l}} = \\mathrm{x}_{k} $$
    %
    % :param sz: size of \\(\\mathrm{y}\\) the output of the :class:`LinOpBroadcast`.
    % :param index: dimensions along which vector will be broadcasted
    %
    % LinOpBroadcast is the Adjoint of the LinOpSum operator
    %
    % All attributes of parent class :class:`LinOp` are inherited.
    %
    % **Example** S=LinOpBroadcast(sz,index)
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %% GUI-Header
    % GUInotation-B-
    % GUIcall-LinOpBroadcast(InputSize,index)-
    % GUIparam-OutputSize-vecInt-[]-Size ofthe output of the LinOpBroadcast (e.g. [512 512])
    % GUIparam-index-vecInt-[]-Dimensions along which vector will be broadcasted (all by default)
    
    %%    Copyright (C) 2018
    %     F. Soulez ferreol.soulez@epfl.ch
    %     Anthony Berdeu (Laboratoire Hubert Curien)
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
    
    %% Constructor
    methods
        function this = LinOpBroadcast(sz,index)
            if nargin == 1
                index = [];
            end
            this.name ='LinOpBroadcast ';
            this.isInvertible=false;
            this.isDifferentiable=true;
            this.sizeout = sz;
            
            this.ndms = length(this.sizeout);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizeout(2) ==1
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
                    this.sizein= [1 1];
                case(this.ndms-1)
                    this.sizein= [this.sizeout(T) 1];
                otherwise
                    this.sizein= this.sizeout(T);
            end
            this.kerdims = this.sizeout;
            this.kerdims(T)=1;
            this.imdims = this.sizeout;
            this.imdims(~T)=1;
            this.norm = prod(this.kerdims); % To be checked
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            % $$\\mathrm{H}^* : \\mathrm{x} \\mapsto \\mathrm{y_{k,l}} =  \\mathrm{x}_{k} \\; \\forall l$$
            y = reshape(repmat(reshape(x,this.imdims),this.kerdims),this.sizeout);
            
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{y_k} = \\sum_l \\mathrm{x}_{k,l} $$
            for n=this.index
                x = sum(x,n);
            end
            y = reshape(x, this.sizein);
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            a = prod(this.kerdims);
            y = x.*a;
        end
        
        function M = makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpSum(this.sizeout, this.index);
        end
        %
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            a = prod(this.kerdims);
            M=LinOpDiag(this.sizein,a);
        end
        
    end
end
