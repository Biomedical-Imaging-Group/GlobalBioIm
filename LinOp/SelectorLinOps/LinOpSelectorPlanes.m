classdef LinOpSelectorPlanes < LinOpSelector
    % Planes Selector linear operator which extracts a subset of planes 
    % from the given volume
    %
    % :param sz:  size of \\(\\mathrm{x}\\) on which the :class:`LinOpSelectorPlane` applies.
    % :param index: dimension along which the planes are extracted      
    % :param pos: vector containing the position of the planes to extract 
    % :param squeezeflag: true to squeeze  singleton dimensions (default: false)
    %
    % All attributes of parent class :class:`LinOpSelector` are inherited. 
    %
    % **Example** S=LinOpSelectorPlanes(sz,index,pos,squeezeflag)
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`,
    
    %%    Copyright (C) 2019
    %     E. Soubies  emmanuel.soubies@irit.fr
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
        index
        pos
        ndms
    end
    
    %% Constructor
    methods
        function this = LinOpSelectorPlanes(sz,index,pos,squeezeflag)
            this.name ='LinOp SelectorPlanes';	
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizein = sz;
            this.ndms = length(this.sizein);
            % Special case for vectors as matlab thought it is matrix ;-(
            if this.sizein(2) ==1
                this.ndms = 1;
            end
            assert(index <= this.ndms,'The index should be a conformable to sz');
            this.index = index;
            assert(isvector(pos) && length(pos)<= this.sizein(this.index) && max(pos)<= this.sizein(this.index),'unvalid sel argument');
            this.pos=pos;
            
            this.sizeout=this.sizein;
            this.sizeout(this.index)=length(this.pos);
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

            this.sel=cell(length(sz),1);
            for ii=1:length(this.sizein)
                if ii==this.index
                    this.sel{ii}=pos;
                else
                    this.sel{ii}=1:this.sizein(ii);
                end
            end
        end
    end
end
