classdef LinOpSum <  LinOp
    % LinOpSum: Summation linear operator which sums the elements of a variable along
    % given directions.
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{y_k} = \\sum_l \\mathrm{x}_{k,l} $$
    %
    % :param sz: size of \\(\\mathrm{x}\\) on which the :class:`LinOpSum` applies.
    % :param index: dimensions along which the sum will be performed (inner sum over l)
    %
    % All attributes of parent class :class:`LinOp` are inherited.
    %
    % **Example** S=LinOpSum(sz,index)
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %% GUI-Header
    % GUInotation-Sum-
    % GUIcall-LinOpSum(InputSize,index)-
    % GUIparam-InputSize-vecInt-[]-Input size of the diagonal operator (e.g. [512 512]).
    % GUIparam-index-vecInt-[]-Dimensions along which the sum will be performed (inner sum over l)
    
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
        index  % index along which dimension are computed the finite differences
        ndms   % number of dimensions of the input
        kerdims % ker dimensions
        imdims % im dimensions
    end
    
    %% Constructor
    methods
        function this = LinOpSum(sz,index)
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            if nargin == 1
                index = [];
            end
            this.name ='LinOpSum ';
            this.isInvertible=false;
            this.isDifferentiable=true;
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
            
            this.norm=sqrt(prod(this.sizein(this.index)));
            
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            for n=this.index
                x = sum(x,n);
            end
            y = reshape(squeeze(x),this.sizeout);
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            % $$\\mathrm{H}^* : \\mathrm{x} \\mapsto \\mathrm{y_{k,l}} =  \\mathrm{x}_{k} \\; \\forall l$$
            y = reshape(repmat(reshape(x,this.imdims),this.kerdims),this.sizein);
        end
        
        function M = makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpBroadcast(this.sizein, this.index);
        end
        %
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            a = prod(this.kerdims);
            y = x.*a;
        end
        function M=makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M= LinOpDiag(this.sizeout,prod(this.kerdims));
        end
        function G = makeComposition_(this, H)
            % Reimplemented from parent class :class:`LinOp`
            
            if isa(H, 'LinOpComposition')
                if isa(H.H2,'LinOpBroadcast') && all(this.kerdims == H.H2.kerdims)
                    if isa(H.H1, 'LinOpConv')
                        idxDiff1=setdiff(this.index,H.H1.Notindex);
                        idxDiff2=setdiff(H.H1.index,this.index);
                        idxUnion=union(this.index,H.H1.Notindex);
                        if isempty(idxDiff2)
                            dd=iSfft(H.H1.mtf,H.H1.Notindex);
                            for n=this.index
                                dd=sum(dd,n);
                            end
                            G=LinOpDiag(H.H2.sizein,squeeze(dd)*prod(this.sizein(idxDiff1)));
                        else
                            newMtf=Sfft(iSfft(H.H1.mtf,H.H1.Notindex),idxUnion);
                            for n=this.index
                                newMtf=sum(newMtf,n);
                            end
                            newMtf=squeeze(newMtf)*prod(this.sizein(idxDiff1));
                            G = LinOpConv(newMtf,H.H1.isReal,idxDiff2);
                        end
                    elseif isa(H.H1,'LinOpDiag')
                        if H.H1.isScaledIdentity
                            a = prod(this.kerdims).*H.H1.diag;
                            G = LinOpDiag(this.sizeout,H.H1.diag);
                        else
                            a = squeeze(sum(H.H1.diag,this.index));
                            G = LinOpDiag(this.sizeout,a);
                        end
                    else
                        G = makeComposition_@LinOp(this, H);
                    end
                else
                    G = makeComposition_@LinOp(this, H);
                end
            elseif isa(H,'LinOpBroadcast') && all(this.kerdims == H.kerdims)
                a = prod(this.kerdims);
                G=LinOpDiag(this.sizeout,a);
            else
                G = makeComposition_@LinOp(this, H);
            end
        end
    end
end
