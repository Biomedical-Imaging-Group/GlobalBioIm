classdef LinOpGrad <  LinOp
    % LinOpGrad:  Gradient linear operator (Finite differences) 
    %
    % :param sz: sizein of the gradient operator
    % :param index: dimension along which the gradient is computed (all by default)
    % :param bc: boundary condition: 'circular' (default), 'zeros', 'mirror'
    % :param res: vector containing the resolution along each dimension (default all 1)
    % :param useRFT: use RFT when defining the :class:`LinOpConv` associated to \\(\\mathrm{H^TH}\\) 
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    % 
    % **Note** When circular boundary conditions are selected, the method
    %          makeHtH (or equivalently the composition ``H'*H``) returns a convolution
    %          linear operator :class:`LinOp`
    %
    % **Example** G = LinOpGrad(sz,index,bc,res)
    %
    % See also :class:`Map`, :class:`LinOp`
    
    %% GUI-Header
    % GUInotation-G-
    % GUIcall-LinOpGrad(InputSize,index,BC,res)-
    % GUIparam-InputSize-vecInt-[]-Input size of the gradient operator (e.g. [512 512]).
    % GUIparam-index-vecInt-[]-Dimensions along which the gradient is computed (all by default)
    % GUIparam-BC-dropDown/circular/zeros/mirror-circular-Boundary condition (default 'circular')
    % GUIparam-res-vecInt-[]-Vector containing the resolution along each dimension (default all 1)
    
    %%    Copyright (C) 2017
    %     F. Soulez ferreol.soulez@epfl.ch, 
    %     E. Soubies emmanuel.soubies@epfl.ch, 
    %     M. McCann michael.mccann@epfl.ch
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
        index;     % index along wich dimension are computed the finite differences
        lgthidx;   % length of INDEX
        ndms;      % number of dimension of the input
        bc;        % boundary condition (default mirror);
        res;       % resolution, vector of lenght ndms
    end
    properties
        useRFT=0;  % use RFT when defining the LinOpConv associated to HtH
    end
    
    %% Constructor
    methods
        function this = LinOpGrad(sz,index,bc,res)
            if nargin == 1, index = [];end
            if nargin<=2 || isempty(bc), bc='circular';end
            if nargin<=3 || isempty(res), res=ones_(size(sz));end
            this.name ='LinOpGrad';
            this.isInvertible=false;
            this.isDifferentiable=true;
            this.res=res;
            this.bc=bc;            
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizein = sz;
            this.ndms = length(this.sizein);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizein(2) ==1
                this.ndms = 1;
            end
            if (~isempty(index))
                assert(isvector(index) && length(index)<= this.ndms && max(index)<= this.ndms,'The index should be a conformable  to sz');
                this.index = index;
            else
                this.index = 1:this.ndms;
            end
            this.lgthidx = length(this.index);
            % size of the output = size of the input x length of the index
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizein(2) ==1
                this.sizeout= [this.sizein(1),this.lgthidx];
            else
                this.sizeout= this.sizein;
            end
            
            if this.lgthidx > 1
                this.sizeout(end+1) = this.lgthidx;
            end
            
            this.norm = 2 * sqrt( sum(1./res.^2) );
            
            validatestring(this.bc, {'mirror', 'circular', 'zeros'});
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
            for diffDimInd = 1:length(this.index)
                diffDim = this.index(diffDimInd);
                
                midElements = allElements;
                midElements{diffDim} = 1:this.sizeout(diffDim)-1;
                y(midElements{:}, diffDimInd)  = diff(x, 1, diffDim)  / this.res(diffDimInd);
                
                lastElement = allElements;
                lastElement{diffDim} = this.sizeout(diffDim);
                switch(this.bc)
                    case('mirror')
                        % y(:,1) = (x([2:end,end])-x)/this.res(1);
                        y(lastElement{:}, diffDimInd)  = 0;
                    case('circular')
                        % y(:,1) = (x([2:end,1])-x)/this.res(1);
                        firstElement = allElements;
                        firstElement{diffDim} = 1;
                        y(lastElement{:}, diffDimInd)  =  (x(firstElement{:}) - x(lastElement{:})) / this.res(diffDimInd);
                    case('zeros')
                        y(lastElement{:}, diffDimInd)  = -x(lastElement{:})  / this.res(diffDimInd);
                end
                
            end
            
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = zeros_(this.sizein);
            allElements = repmat({':'}, 1, this.ndms);
            for diffDimInd = 1:length(this.index)
                diffDim = this.index(diffDimInd);
                midElements = allElements;
                midElements{diffDim} = 2:this.sizein(diffDim)-1;
                leftElements = allElements;
                leftElements{diffDim} = 1:this.sizein(diffDim)-1;
                y(midElements{:})  = y(midElements{:}) + diff(-x(leftElements{:}, diffDimInd), 1, diffDim) / this.res(diffDimInd);
                % handle boundary conditions
                lastElement = allElements;
                lastElement{diffDim} = this.sizeout(diffDim);
                firstElement = allElements;
                firstElement{diffDim} = 1;
                secondToLastElement = allElements;
                secondToLastElement{diffDim} = this.sizeout(diffDim) - 1;
                switch(this.bc)
                    case('mirror')
                        %y= [[-x(1) ; (-x(2:end-1)+x(1:end-2))] ; x(end-1)]/this.res(1);
                        y(firstElement{:})  =  y(firstElement{:}) + -x(firstElement{:}, diffDimInd)/this.res(diffDimInd);
                        y(lastElement{:})   = y(lastElement{:}) + x(secondToLastElement{:}, diffDimInd)/this.res(diffDimInd);
                    case('circular')
                        %y= [[x(end)-x(1) ; (-x(2:end-1)+x(1:end-2))] ; (x(end-1)-x(end))]/this.res(1);
                        y(firstElement{:})  =  y(firstElement{:}) ...
                            + (x(lastElement{:}, diffDimInd)-x(firstElement{:}, diffDimInd))/this.res(diffDimInd);
                        y(lastElement{:}) = y(lastElement{:}) ...
                            + (x(secondToLastElement{:}, diffDimInd) -  x(lastElement{:}, diffDimInd))/this.res(diffDimInd);
                    case('zeros')
                        %y= [[-x(1) ; (-x(2:end-1)+x(1:end-2))] ; x(end-1)-x(end)]/this.res(1);
                        y(firstElement{:})  =  y(firstElement{:}) ...
                            + (-x(firstElement{:}, diffDimInd))/this.res(diffDimInd);
                        y(lastElement{:}) = y(lastElement{:}) ...
                            + (x(secondToLastElement{:}, diffDimInd) -  x(lastElement{:}, diffDimInd))/this.res(diffDimInd);
                end
            end
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = zeros_(this.sizein);
            allElements = repmat({':'}, 1, this.ndms);
            for diffDimInd = 1:length(this.index)
                diffDim = this.index(diffDimInd);
                midElements = allElements;
                midElements{diffDim} = 2:this.sizein(diffDim)-1;
                y(midElements{:}) = y(midElements{:}) + -diff(x, 2, diffDim)/this.res(diffDimInd)^2;
                
                lastElement = allElements;
                lastElement{diffDim} = this.sizeout(diffDim);
                
                firstElement = allElements;
                firstElement{diffDim} = 1;
                
                secondElement = allElements;
                secondElement{diffDim} = 2;
                
                secondToLastElement = allElements;
                secondToLastElement{diffDim} = this.sizeout(diffDim) - 1;
                switch(this.bc)
                    case('mirror')
                        % y = (2*x - x([1,1:end-1]) - x([2:end,end]))/this.res(1)^2;
                        y(firstElement{:}) = y(firstElement{:}) ...
                            + ( x(firstElement{:}) - x(secondElement{:}) ) / this.res(diffDimInd)^2;
                        y(lastElement{:}) = y(lastElement{:}) ...
                            + (x(lastElement{:}) - x(secondToLastElement{:})) / this.res(diffDimInd)^2;
                    case('circular')
                        % y = (2*x - x([end,1:end-1]) - x([2:end,1]))/this.res(1)^2;
                        y(firstElement{:}) = y(firstElement{:}) ...
                            + ( 2*x(firstElement{:}) - x(lastElement{:}) - x(secondElement{:}) ) / this.res(diffDimInd)^2;
                        y(lastElement{:}) = y(lastElement{:}) ...
                            + (2*x(lastElement{:}) - x(secondToLastElement{:}) - x(firstElement{:}) ) / this.res(diffDimInd)^2;
                    case('zeros')
                        % y(1:end-1) = (2*x(1:end-1) - x([1,1:end-2]) - x(2:end))/this.res(1)^2;
                        % y(end)=(2*x(end)-x(end-1))/this.res(1)^2;
                        y(firstElement{:}) = y(firstElement{:}) ...
                            + ( x(firstElement{:}) - x(secondElement{:}) ) / this.res(diffDimInd)^2;
                        y(lastElement{:}) = y(lastElement{:}) ...
                            + (2*x(lastElement{:}) - x(secondToLastElement{:}) ) / this.res(diffDimInd)^2;
                end
            end
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if strcmp(this.bc,'circular')
                fHtH=zeros_(this.sizein)  ;
                idxAll = repmat({1}, 1, max(this.ndms,2));  % The  max(this.ndms,2) is to deal with the vectorial case (i.e. this.ndms==1)
                rep=this.sizein;
                sel=idxAll;
                for ii = this.index
                    idx=idxAll;idx{ii}=':';
                    dd=idxAll;dd{ii}=this.sizein(ii);
                    fHtH(idx{:})=fHtH(idx{:})+reshape([2, -1, zeros(1,this.sizein(ii)-3),  -1],dd{:})/this.res(ii)^2;
                    rep(ii)=1;sel{ii}=':';
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
