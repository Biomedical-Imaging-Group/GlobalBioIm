classdef LinOpSumPatches <  LinOp
    % LinOpSumPatches: Linear operator which sums defined patches of a
    % variable.
    %
    % :param sz: size of \\(\\mathrm{x}\\) on which the :class:`LinOpSum` applies.
    % :param szPatch: array containing the patch size in each direction
    % 
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** S=LinOpSumPatches(sz,szPatches)
    %
    % See also :class:`LinOp`, :class:`Map`
    
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
        szPatch    % array containing the patch size in each direction
        sel;
    end
    
    %% Constructor
    methods
        function this = LinOpSumPatches(sz,szPatch)
            this.name ='LinOpSumPatches';
            assert(cmpSize(size(sz),size(szPatch)),'Parameters sz and szPatches must have the same size');
            assert(~any(mod(sz,szPatch)),'Sizes in sz must be multiples of patch sizes in szPatches');  
            this.sizein=sz;
            this.szPatch=szPatch;
            this.sizeout=this.szPatch;           
            for n=1:length(this.sizein)
                this.sel{n}=this.szPatch(n)*ones(1,this.sizein(n)/this.szPatch(n));
            end
		end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.  
            y=zeros_(this.szPatch);
            tmp=mat2cell(x,this.sel{:});
            for n=1:numel(tmp)
                y=y+tmp{n};
            end
        end		
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=repmat(x,this.sizein./this.szPatch);
        end		
    end
end