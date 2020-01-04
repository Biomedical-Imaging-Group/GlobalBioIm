classdef LinOpHess <  LinOp
    % LinOpHess: Hessian linear operator (finite differences) 
    %
    % :param sz: sizein of the gradient operator
    % :param bc: boundary condition: 'circular' (default), 'mirror'
    % :param useRFT: use RFT when defining the :class:`LinOpConv` associated to \\(\\mathrm{H^TH}\\)     
    % :param index: dimension along which the hessian is computed (all by default)
    %
    % The output size is [sz,lidx(lidx+1)/2] where lidx is the length of
    % the index vector. The last dimension of length lidx(lidx+1)/2
    % represents the coefficients of upper-right part of the symmetric hessian 
    % matrix ordered from top to bottom and left to right. For instance:
    %
    % - in 3D with index=[1 2 3], [d^2F/dxx;d^2F/dxy;d^2F/dxz;d^2F/dyy;d^2F/dyz;d^2F/dzz]
    %
    % - in 3D with index=[2 3], [d^2F/dyy;d^2F/dyz;d^2F/dzz]
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note** When circular boundary conditions are selected, the method
    % makeHtH (or equivalently the composition ``H'*H``) returns a convolution
    % linear operator :class:`LinOp`
    %
    % **Example** H = LinOpHess(sz,bc,index)
    %
    % See also :class:`Map`, :class:`LinOp`
    
    %% GUI-Header
    % GUInotation-H-
    % GUIcall-LinOpHess(InputSize,BC,index)-
    % GUIparam-InputSize-vecInt-[]-Input size of the gradient operator (e.g. [512 512]).
    % GUIparam-index-vecInt-[]-dimension along which the gradient is computed (all by default)
    % GUIparam-BC-dropDown/circular/mirror-circular-Boundary condition (default 'circular')
    
	%%    Copyright (C) 2017 
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
	  ndms;      % number of dimension of the input
      bc;        % boundary condition (default mirror);
      index;     % index along wich dimension are computed the hessian    
      lgthidx;   % length of INDEX
    end
    properties
        useRFT=0;  % use RFT when defining the LinOpConv associated to HtH
    end
    
    %% Constructor
    methods
        function this = LinOpHess(sz,bc,index)
            if nargin < 3, index=[];end
            if nargin < 2, index=[]; bc=[];end
            if isempty(bc), bc='circular'; end
            if isempty(index), index=1:length(sz); end
            this.name='LinOpHess';
            this.isInvertible=false;
            this.isDifferentiable=true;
            this.sizein=sz;
            this.ndms = length(this.sizein);
            this.bc=bc;
            this.index=index;
            this.lgthidx = length(this.index);
            this.sizeout= this.sizein;
            if this.lgthidx > 1
                this.sizeout(end+1) = (this.lgthidx+1)*this.lgthidx/2;
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
    % - applyHtH_(this,x)
    % - makeHtH_(this)
    methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = zeros_(this.sizeout);
            allElements = repmat({':'}, 1, this.ndms);
            idx=1;
            for DimInd = 1:length(this.index)
                Dim = this.index(DimInd);
                ll = allElements;rr = allElements;
                % Second derivative for DimInd
                switch(this.bc)
                    case('circular')
                        rr{Dim} = [3:this.sizeout(Dim),1,2];
                        ll{Dim} = [2:this.sizeout(Dim),1];
                        y(allElements{:}, idx)=x(rr{:}) -2*x(ll{:}) +x;
                    case('mirror')
                        rr{Dim} = [3:this.sizeout(Dim),this.sizeout(Dim),this.sizeout(Dim)-1];
                        ll{Dim} = [2:this.sizeout(Dim),this.sizeout(Dim)];
                        y(allElements{:}, idx)=x(rr{:}) -2*x(ll{:}) +x;
                end
                idx=idx+1;
                % Cross derivatives for DimInd and remaining dimensions
                for DimInd2=DimInd+1:length(this.index)
                    Dim2 = this.index(DimInd2);
                    DimCrossElem=ll;
                    switch(this.bc)
                        case('circular')
                            DimCrossElem{Dim2} = [2:this.sizeout(Dim2),1];
                            ll2=allElements; ll2{Dim2} = [2:this.sizeout(Dim2),1];
                            y(allElements{:}, idx)=x(DimCrossElem{:}) - x(ll{:}) - x(ll2{:}) +x;
                        case('mirror')
                            DimCrossElem{Dim2} = [2:this.sizeout(Dim2),this.sizeout(Dim2)];
                            ll2=allElements; ll2{Dim2} = [2:this.sizeout(Dim2),this.sizeout(Dim2)];
                            y(allElements{:}, idx)=x(DimCrossElem{:}) - x(ll{:}) - x(ll2{:}) +x;
                    end
                    idx=idx+1;
                end
            end
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`            
            y = zeros_(this.sizein);
            allElements = repmat({':'}, 1, this.ndms);
            idx=1;
            for DimInd = 1:length(this.index)
                Dim = this.index(DimInd);
                switch(this.bc)
                    case('circular')
                        elem1=allElements;elem1{Dim}=[this.sizein(Dim)-1,this.sizein(Dim),1:this.sizein(Dim)-2];
                        elem2=allElements;elem2{Dim}=[this.sizein(Dim),1:this.sizein(Dim)-1];
                        y=y + x(elem1{:},idx) -2*x(elem2{:},idx) + x(allElements{:},idx);
                    case('mirror')
                        elem1=allElements;elem2=allElements;elem3=allElements;elem4=allElements;
                        elem1{Dim}=1; elem2{Dim}=2;
                        y(elem1{:})=y(elem1{:}) + x(elem1{:},idx);
                        y(elem2{:})=y(elem2{:}) + x(elem2{:},idx) - 2*x(elem1{:},idx);
                        elem1{Dim}=3:this.sizein(Dim)-2;elem2{Dim}=1:this.sizein(Dim)-4;
                        elem3{Dim}=2:this.sizein(Dim)-3;
                        y(elem1{:})=y(elem1{:}) + x(elem2{:},idx) -2*x(elem3{:},idx) + x(elem1{:},idx);
                        elem1{Dim}=this.sizein(Dim)-1;elem2{Dim}=this.sizein(Dim)-3;
                        elem3{Dim}=this.sizein(Dim)-2;elem4{Dim}=this.sizein(Dim);
                        y(elem1{:})=y(elem1{:}) + x(elem2{:},idx) -2*x(elem3{:},idx) + x(elem1{:},idx) + x(elem4{:},idx);
                        y(elem4{:})=y(elem4{:}) + x(elem3{:},idx) - x(elem1{:},idx) - x(elem4{:},idx);
                end
                idx=idx+1;
                for DimInd2=DimInd+1:length(this.index)
                    Dim2 = this.index(DimInd2);
                    switch(this.bc)
                        case('circular')
                            DimCrossElem=elem2;
                            DimCrossElem{Dim2} =[this.sizein(Dim2),1:this.sizein(Dim2)-1];
                            elem22=allElements;elem22{Dim2} =[this.sizein(Dim2),1:this.sizein(Dim2)-1];
                            y= y + x(DimCrossElem{:},idx) - x(elem22{:},idx)  - x(elem2{:},idx) + x(allElements{:},idx);
                        case('mirror')
                            elem1=allElements;elem2=allElements;elem3=allElements;elem4=allElements;
                            elem1{Dim}=2:this.sizein(Dim)-1;elem1{Dim2}=2:this.sizein(Dim2)-1;
                            elem2{Dim}=1:this.sizein(Dim)-2;elem2{Dim2}=1:this.sizein(Dim2)-2;
                            elem3{Dim}=1:this.sizein(Dim)-2;elem3{Dim2}=2:this.sizein(Dim2)-1;
                            elem4{Dim}=2:this.sizein(Dim)-1;elem4{Dim2}=1:this.sizein(Dim2)-2;
                            y(elem1{:})=y(elem1{:}) + x(elem2{:},idx) - x(elem3{:},idx) - x(elem4{:},idx) +x(elem1{:},idx);
                            elem1{Dim2}=1;elem2{Dim2}=1;
                            y(elem1{:})=y(elem1{:}) -x(elem2{:},idx) + x(elem1{:},idx);
                            elem3{Dim}=1;elem4{Dim}=1;
                            y(elem3{:})=y(elem3{:}) -x(elem4{:},idx) + x(elem3{:},idx);
                            elem1{Dim2}=this.sizein(Dim2);
                            elem2{Dim}=2:this.sizein(Dim)-1;elem2{Dim2}=this.sizein(Dim2)-1;
                            elem3{Dim}=1:this.sizein(Dim)-2;elem3{Dim2}=this.sizein(Dim2)-1;
                            y(elem1{:})=y(elem1{:}) -x(elem2{:},idx) + x(elem3{:},idx);
                            elem1{Dim}=this.sizein(Dim);elem1{Dim2}=2:this.sizein(Dim2)-1;
                            elem2{Dim}=this.sizein(Dim)-1;elem2{Dim2}=2:this.sizein(Dim2)-1;
                            elem3{Dim}=this.sizein(Dim)-1;elem3{Dim2}=1:this.sizein(Dim2)-2;
                            y(elem1{:})=y(elem1{:}) -x(elem2{:},idx) + x(elem3{:},idx);
                            elem1{Dim}=1;elem1{Dim2}=1;
                            y(elem1{:})=y(elem1{:})+x(elem1{:},idx);
                            elem1{Dim}=this.sizein(Dim);elem1{Dim2}=this.sizein(Dim2);
                            elem2{Dim}=this.sizein(Dim)-1;elem2{Dim2}=this.sizein(Dim2)-1;
                            y(elem1{:})=y(elem1{:})+x(elem2{:},idx);
                            elem1{Dim}=1;elem1{Dim2}=this.sizein(Dim2);
                            elem2{Dim}=1;elem2{Dim2}=this.sizein(Dim2)-1;
                            y(elem1{:})=y(elem1{:})-x(elem2{:},idx);
                            elem1{Dim}=this.sizein(Dim);elem1{Dim2}=1;
                            elem2{Dim}=this.sizein(Dim)-1;elem2{Dim2}=1;
                            y(elem1{:})=y(elem1{:})-x(elem2{:},idx);
                    end
                    idx=idx+1;
                end
            end 
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if strcmp(this.bc,'circular')
                fHtH=zeros_(this.sizein)  ;
                idxAll = repmat({1}, 1, this.ndms);
                rep=this.sizein;
                sel=idxAll;
                for ii = 1:length(this.index)
                    Dim = this.index(ii);
                    idx=idxAll;idx{Dim}=':';
                    dd=idxAll;dd{Dim}=this.sizein(Dim);
                    fHtH(idx{:})=fHtH(idx{:})+reshape([6, -4, 1, zeros(1,this.sizein(Dim)-5), 1, -4],dd{:});
                    for jj=ii+1:length(this.index)
                        Dim2 = this.index(jj);
                        idx=idxAll;idx{Dim}=[1 2];idx{Dim2}=[1 2];
                        dd=idxAll;dd{Dim}=2;dd{Dim2}=2;
                        fHtH(idx{:})=fHtH(idx{:})+reshape([4, -2;-2 1],dd{:});
                        idx=idxAll;idx{Dim}=[1 2];idx{Dim2}=this.sizein(Dim2);
                        dd=idxAll;dd{Dim}=2;dd{Dim2}=1;
                        fHtH(idx{:})=fHtH(idx{:})+reshape([-2 1],dd{:});
                        idx=idxAll;idx{Dim}=this.sizein(Dim);idx{Dim2}=[1 2];
                        dd=idxAll;dd{Dim}=1;dd{Dim2}=2;
                        fHtH(idx{:})=fHtH(idx{:})+reshape([-2 1],dd{:});
                        idx=idxAll;idx{Dim}=this.sizein(Dim);idx{Dim2}=this.sizein(Dim2);
                        dd=idxAll;dd{Dim}=1;dd{Dim2}=1;
                        fHtH(idx{:})=fHtH(idx{:})+1;
                    end
                    rep(Dim)=1;sel{Dim}=':';
                end
                fHtH=repmat(fHtH(sel{:}) ,rep);
                if this.useRFT
                    M=LinOpConv('PSF',fHtH,1,this.index,'useRFT');
                else
                    M=LinOpConv('PSF',fHtH,1,this.index);
                end
            else
                M=makeHtH_@LinOp(this);
            end
        end
    end
end

