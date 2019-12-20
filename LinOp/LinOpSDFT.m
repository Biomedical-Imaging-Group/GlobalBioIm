classdef LinOpSDFT <  LinOp
    % LinOpSDFT : Sliced Discrete Fourier operator
    %
    % :param sz: sizein of the operator.
    % :param unitary: boolean true when normalized DFT (default false)
    % :param pad: padding size (see the doc of fftn function).
    % :param index: index along wich dimension are computed the FFT (default all)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** SDFT=LinOpSDFT(sz,index,unitary, pad)
    %
    % See also :class:`LinOp`, :class:`Map`, fftn, ifftn, Sfft, iSfft
    
    %% GUI-Header
    % GUInotation-S-
    % GUIcall-LinOpSelectorPatch(InputSize,index,unitary,pad)-
    % GUIparam-InputSize-vecInt-[]-Input size of the selector operator
    % GUIparam-index-vecInt-[]-Dimensions along which the DFT is computed (default all)
    % GUIparam-unitary-boolean-0-boolean true when normalized DFT (default false)
    % GUIparam-pad-vecInt-[]-padding size (see the doc of fftn function).

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
        index % index along wich dimension are computed the FFT
        Notindex% ~index
        N % Number of element
        ndms % number of dimensions
        unitary
    end
    
    %% Constructor
    methods
        function this = LinOpSDFT(sz,index,unitary, pad)
            if nargin < 2 || isempty(index), index= 1:length(sz); end;
            if nargin < 3 || isempty(unitary), unitary=false; end;
            if nargin < 4 || isempty(pad), pad=sz; end
            assert(issize(index),'The index should be a conformable  to sz');
            assert(islogical(unitary), 'UNITARY should be logical');  
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            assert(issize(pad) || isempty(pad),'pad  should be a conformable  to a size ');
            this.index = index;
            this.unitary= unitary;
            this.sizeout=pad;
            this.sizein=sz;
            this.name ='LinOpSDFT';
            if all(this.sizein==this.sizeout)
                this.isInvertible=true; % i.e. there is no padding
            end
            this.isDifferentiable=true;            
            this.ndms = length(this.sizein);
            % Special case for vectors 
            if (this.ndms==2) && (this.sizein(2) ==1 || this.sizein(1) ==1)
                this.ndms = 1;
            end           
            if (~isempty(this.index))
                dim = 1:this.ndms;
                Iidx = true(this.ndms,1);
                Iidx(this.index) = 0;
                this.Notindex = dim(Iidx);
            else
                this.index = 1:this.ndms;
                this.Notindex = [];
            end
            assert(isequal(this.sizein(this.Notindex),pad(this.Notindex)),'padding values associated to NotIndex elements should be equal to those of sizein');
            this.N= prod(this.sizeout(this.index));    
            if this.unitary
                this.norm=1;
            else
                this.norm=sqrt(this.N);
            end               
        end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                y = 1./sqrt(this.N) * Sfft(x,this.Notindex,this.sizeout);
            else
                y =  Sfft(x,this.Notindex,this.sizeout);
            end
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                y = sqrt(this.N) * iSfft(x,this.Notindex,this.sizein);
            else
                y = this.N * iSfft(x,this.Notindex,this.sizein);
            end
        end
        function y = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                if this.unitary
                    y =  iSfft(x*sqrt(this.N), this.Notindex);
                else
                    y = iSfft(x, this.Notindex);
                end
            else
                y=applyInverse_@LinOp(this,x);
            end
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                y = x;
            else
                y =this.N *x;
            end
        end
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if all(this.sizein==this.sizeout)
                y=this.applyHtH_(x);
            else
                y=applyHHt_@LinOp(this,x);
            end
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`. 
            if this.isInvertible
                if this.unitary
                    y = 1/sqrt(this.N) * Sfft(x, this.Notindex);
                else
                    y = 1/this.N * Sfft(x, this.Notindex);
                end
            else
                y=applyAdjointInverse_@LinOp(this,x);
            end
        end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if all(this.sizein==this.sizeout)
                if this.unitary
                    M=LinOpDiag(this.sizein);
                else
                    M=LinOpDiag(this.sizein,this.N);
                end
            else
                M=makeHHt_@LinOp(this);
            end
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                M=LinOpDiag(this.sizein);
            else
                M=LinOpDiag(this.sizein,this.N );
            end
        end
        function M = makeInversion_(this)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                if this.unitary
                    M = this.makeAdjoint;
                else
                    M = 1./this.N * this.makeAdjoint;
                end
            else
                M=makeInversion_@LinOp(this);
            end
        end
    end
end

