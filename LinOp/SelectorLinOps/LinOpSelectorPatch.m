classdef LinOpSelectorPatch < LinOpSelector
    % Patch Selector linear operator which extracts a patch from the given
    % vector
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{x}_{[i_{min}:i_{max}]}$$
    % where \\( i_{min} \\) and \\( i_{max} \\) are indexes corresponding
    % to the first and last elements of the patch.
    %
    % :param sz:  size of \\(\\mathrm{x}\\) on which the :class:`LinOpSelectorPatch` applies.
    % :param idxmin: array containing the first kept index in each direction         
    % :param idxmax: array containing the last kept index in each direction  
    % :param squeezeflag: true to squeeze  singleton dimensions (default: false)
    %
    % All attributes of parent class :class:`LinOpSelector` are inherited. 
    %
    % **Example** S=LinOpSelectorPatch(sz,idxmin,idxmax)
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`,
    % :class:`LinOpDownsampling`.
    
    %% GUI-Header
    % GUInotation-S-
    % GUIcall-LinOpSelectorPatch(InputSize,idxmin,idxmax,1)-
    % GUIparam-InputSize-vecInt-[]-Input size of the selector operator
    % GUIparam-idxmin-vecInt-[]-Array containing the first kept index in each direction
    % GUIparam-idxmax-vecInt-[]-Array containing the last kept index in each direction
    
    %%    Copyright (C) 2017 
    %     E. Soubies  emmanuel.soubies@epfl.ch
    %     F. Soulez 
    %     A. Berdeu
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
        function this = LinOpSelectorPatch(sz,idxmin,idxmax,squeezeflag)
            this.name ='LinOp SelectorPatch';	
            assert(cmpSize(size(sz),size(idxmin)) && cmpSize(size(sz),size(idxmax)),'Parameters sz idxmin and idxmax must have the same size');
            assert(~any(idxmin<=0) && ~any(idxmin>sz),'idxmin out of sz range');
            assert(~any(idxmax<=0) && ~any(idxmax>sz),'idxmax out of sz range');
            assert(~any(idxmin>idxmax),'idxmin must be smaller than idxmax in each dimension');
            this.idxmin=idxmin;
            this.idxmax=idxmax;
            this.sizeout=idxmax-idxmin+1;
            % flag_squeeze
            if (nargin>3)&& ~isscalar(squeezeflag )
                error('squeezeflag must be a single boolean');
            end
            if (nargin>3) && squeezeflag
                this.sizeout = this.sizeout(this.sizeout~=1);
            end
            if numel(this.sizeout)==1
                this.sizeout = [this.sizeout, 1];
            end
            this.sizein=sz;
            this.sel=cell(length(sz),1);
            for ii=1:length(sz)
                this.sel{ii}=this.idxmin(ii):this.idxmax(ii);
            end
        end
    end
end
