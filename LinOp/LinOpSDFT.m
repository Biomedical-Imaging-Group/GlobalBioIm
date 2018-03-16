classdef LinOpSDFT <  LinOp
    % LinOpSDFT : Sliced Discrete Fourier operator
    %
    % :param sz: sizein of the operator.
    % :param index: index along wich dimension are computed the FFT (default all)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** SDFT=LinOpSDFT(sz,index)
    %
    % See also :class:`LinOp`, :class:`Map`, :class:`LinOpDFT`, fftn,
    % ifftn, Sfft, iSfft
    
    % FIXME : Should better be merged with DFT

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
    end
    
    %% Constructor
    methods
        function this = LinOpSDFT(sz,index)
            if nargin < 2 || (isempty(index))
                this.index = 1:length(sz);
            else
                assert(issize(index),'The index should be a conformable  to sz');
                this.index = index;
            end           
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizeout=sz;
            this.sizein=sz;

            this.name ='LinOpSDFT';
            this.isInvertible=true;
            this.isDifferentiable=true;
            
            this.ndms = length(this.sizein);
            % Special case for vectors as matlab thought it is matrix ;-(
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
            this.N= prod(this.sizein(this.index));
        end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =  Sfft(x, this.Notindex);
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =   this.N * iSfft(x,this.Notindex);
        end
        function y = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =  iSfft(x, this.Notindex);
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = this.N * x;
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =   1/this.N * Sfft(x, this.Notindex);
        end
    end
end

