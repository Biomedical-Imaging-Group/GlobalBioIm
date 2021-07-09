classdef LinOpPad <  LinOpSelector
    % Selector linear operator which pad the input vector with zeros 
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{x}_{\\mathrm{padded}} $$
    % where \\(\\mathrm{reg}\\) is a region of the input vector.
    %
    % It is the adjoint of LinOpCrop
    %
    % :param sizein:  size of \\(\\mathrm{x}\\) on which the :class:`LinOpPad` applies.
    % :param padsize: number of zeros before and after the input vector: 
    %                 padsize(1,n): number of zeros before the first element along dimension n
    %                 padsize(2,n): number of zeros after the last element along dimension n  
    %
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    %
    % **Example** S=LinOpPad([128,256],[ [128 256 ]' [128,128 ]'])
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`, :class:`LinOpCrop`
    % :class:`LinOpDownsample`.
    
    %%    Copyright (C) 2020
    %     F. Soulez ferreol.soulez@univ-lyon1.fr 
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
        padsize;         % cell containing a boolean array (the cell serve for consistency with derived classes)
        str;
    end
    
    %% Constructor
    methods
        function this = LinOpPad(sizein,padsize)
            this.name ='LinOpPad';
            assert(issize(sizein),'sizein must be a size');
            this.sizein=sizein;
            assert(size(padsize,1)==2,'padsize first dim must be 2');
            str='sel(';
            for n=1:numel(sizein)
                if n<=size(padsize,2)
                    assert( (padsize(2,n)>=0) && (padsize(1,n) >=0),['padsize(,',n,') must be positive.']);
                else
                    padsize(1,n)=1;
                    padsize(2,n)=sizein(n);                    
                end
                this.sizeout(n) = padsize(2,n) + padsize(1,n)+sizein(n);
                str = [str,num2str(padsize(1,n)+1),':',num2str(padsize(1,n)+sizein(n))];
                if n < numel(sizein)
                    str = [str,','];
                end
            end
            this.str= [str,')=true;'];
            sel = logical(zeros_(this.sizeout));
            this.padsize=padsize;
            eval(this.str);
            this.sel = {sel};
            this.norm = 1;	
            this.isInvertible=false;
        end    
    end
      %% Core Methods containing implementations (Protected)
    methods (Access = protected)	
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.           
            y =x(this.sel{:});
            y = reshape(y, this.sizein);
        end        
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.
            y = zeros_(this.sizeout);
            y(this.sel{:}) = x;
        end
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.
            y = zeros_(this.sizeout);
            y(this.sel{:}) = x(this.sel{:});
        end
        function y = applyHtH_(~,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            y = x;
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOpSelector`.
            M=LinOpIdentity(this.sizein);
        end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOpSelector`.
            w=zeros_(this.sizeout);w(this.sel{:})=1;
            M=LinOpDiag([],w);
        end
        
        function M = makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOp`.
            reg = this.padsize;
            reg(2,:) = reg(1,:)+this.sizein;
            reg(1,:) = reg(1,:)+1;
            M=LinOpCrop(this.sizeout, reg);
        end
    end
  
end
