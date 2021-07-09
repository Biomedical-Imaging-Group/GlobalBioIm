classdef LinOpCrop <  LinOpSelector
    % Selector linear operator which crop a specific region of a vector
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{x}_{\\mathrm{reg}} $$
    % where \\(\\mathrm{reg}\\) is a region of the input vector.
    %
    % :param sizein:  size of \\(\\mathrm{x}\\) on which the :class:`LinOpCrop` applies.
    % :param region:  selected region = [[dim1min dim2min ... dimNmin];
    %                                          [dim1max dim2max ... dimNmax]]
    %                  indices can be negative
    % All attributes of parent class :class:`LinOp` are inherited.
    %
    %
    % **Example** S=LinOpCrop(sizein,reg)
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`,:class:`LinOpPad`
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
        region;         % cell containing a boolean array (the cell serve for consistency with derived classes)
        str;
    end
    
    %% Constructor
    methods
        function this = LinOpCrop(sizein,region)
            this.name ='LinOpCrop';
            assert(issize(sizein),'sizein must be a size');
            this.sizein=sizein;
            assert(size(region,1)==2,'region first dim must be 2');
            this.sizeout = sizein;
            sel = logical(zeros_(sizein));
            str='sel(';
            for n=1:numel(sizein)
                if n<=size(region,2)
                    if region(1,n)<0
                        region(1,n) = sizein(n) + region(1,n);
                    end
                    if region(2,n)<0
                        region(2,n) = sizein(n) + region(2,n);
                    end
                    assert( region(2,n)-region(1,n) >=0,['region(2,',n,') must be higher than region(1,',n,').']);
                else
                    region(1,n)=1;
                    region(2,n)=sizein(n);
                end
                this.sizeout(n) = region(2,n) - region(1,n)+1;
                str = [str,num2str(region(1,n)),':',num2str(region(2,n))];
                if n < numel(sizein)
                    str = [str,','];
                end
            end
            this.str= [str,')=true;'];
            this.region=region;
            eval(this.str);
            this.sel = {sel};
            this.norm = 1;
            this.isInvertible=false;
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        
        function M = makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOpSelector`.
            padsize = this.region;
            padsize(2,:) = this.sizein - padsize(2,:);
            padsize(1,:) = padsize(1,:)-1;
            M=LinOpPad(this.sizeout, padsize);
        end
    end
end
