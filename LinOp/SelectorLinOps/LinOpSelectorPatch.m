classdef LinOpSelectorPatch < LinOpSelector
    % Patch Selector linear operator which extracts a patch from the given
    % vector
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{x}_{[i_{min}:i_{max}]}$$
    % where \\( i_{min} \\) and \\( i_{max} \\) are indexes corresponding
    % to the first and last elements of the patch.
    %
    % :param sz:  size of \\(\\mathrm{x}\\) on which the :class:`LinOpDownsample` applies.
    % :param idxmin: array containing the first kept index in each direction         
    % :param idxmax: array containing the last kept index in each direction  
    %
    % All attributes of parent class :class:`LinOpSelector` are inherited. 
    %
    % **Example** S=LinOpSelectorPatch(sz,idxmin,idxmax)
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`,
    % :class:`LinOpDownsampling`.
    
    %%    Copyright (C) 2017 
    %     E. Soubies  emmanuel.soubies@epfl.ch
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
        idxmin
        idxmax
    end
    
    %% Constructor
    methods
        function this = LinOpSelectorPatch(sz,idxmin,idxmax)
            this.name ='LinOp SelectorPatch';	
            assert(isequal(size(sz),size(idxmin)) && isequal(size(sz),size(idxmax)),'Parameters sz idxmin and idxmax must have the same size');
            assert(~any(idxmin<=0) && ~any(idxmin>sz),'idxmin out of sz range');
            assert(~any(idxmax<=0) && ~any(idxmax>sz),'idxmax out of sz range');
            assert(~any(idxmin>idxmax),'idxmin must be smaller than idxmax in each dimension');
            this.idxmin=idxmin;
            this.idxmax=idxmax;
            this.sizeout=idxmax-idxmin+1;
            if this.sizeout(end)==1, this.sizeout=this.sizeout(1:end-1);end;
            this.sizein=sz;
            this.sel=cell(length(sz),1);
            for ii=1:length(sz)
                this.sel{ii}=this.idxmin(ii):this.idxmax(ii);
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)		
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.           
            assert(isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y =x(this.sel{:});
        end        
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            y = zeros_(this.sizein);
            y(this.sel{:}) = x;
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = zeros_(this.sizein);
            y(this.sel{:}) = x(this.sel{:});            
        end
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            y = x;
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOpSelector`. 
            w=zeros_(this.sizein);w(this.sel{:})=1;
            M=LinOpDiag([],w);
        end
    end
end