classdef LinOpDFT <  LinOp
    % LinOpDFT : Discrete Fourier operator
    %
    % Compute the dicrete fourier transfom padded by PAD (optionnal) such
    % that (equivalent to y = fftn(x,pad))
    %
    % :param sz: sizein of the operator.
    % :param pad: padding size (see the doc of fftn function).
    % :param unitary: boolean true when normalized DFT (default false)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** DFT=LinOpDFT(sz, pad)
    %
    % See also :class:`LinOp`, :class:`Map`, :class:`LinOpSDFT`, fftn,
    % ifftn
    
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
        N
        pad  % size of of the fourier (default all)
        real % 1 if real complex Fourier transform (default 0).  FIXME: NOT IMPLEMENTED
        unitary = false
    end
    
    %% Constructor
    methods
        function this = LinOpDFT(sz, pad)        
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.name ='LinOpDFT';
            this.isInvertible=true;
            this.isDifferentiable=true;           
            this.sizein = sz;
            this.sizeout=sz;
            this.N=prod(sz);
            if nargin>1
                if issize(pad)
                    this.pad= pad  ;
                else
                    error('PAD should be conformable to size() output');
                end
            end           
		end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                y = 1./sqrt(this.N) * fftn(x,this.pad);
            else
                y =  fftn(x,this.pad);              
            end
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                y = sqrt(this.N) * ifftn(x,this.pad);
            else
                y = this.N * ifftn(x,this.pad);
            end
        end
        function y = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                y = ifftn(x*sqrt(this.N),this.pad);
            else
                y = ifftn(x,this.pad);
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
        function y = applyAdjointInverse_(this,x) 
            % Reimplemented from parent class :class:`LinOp`.
            if this.unitary
                y = 1/sqrt(this.N) * fftn(x,this.pad);
            else
                y = 1/this.N * fftn(x,this.pad);
            end
        end
    end
end
