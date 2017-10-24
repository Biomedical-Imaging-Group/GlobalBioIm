classdef LinOpDownsample < LinOpSelector
    % Downsampling linear operator 
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{y}_k =\\mathrm{x}_{(k-1)d+1} $$
    % where \\(d \\) is the downsampling factor.
    %
    % :param sz:  size of \\(\\mathrm{x}\\) on which the :class:`LinOpDownsample` applies.
    % :param df: array containing the downsampling factor in each direction
    %            (must divide the corresponding entry of sz)
    %
    % All attributes of the parent class :class:`LinOpSelector` 
    % are inherited. 
    %
    % **Example** D=LinOpDownsample(sz,df)
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`,
    % :class:`LinOpSelectorPatch`.
    
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
        df; % downsampling factors in each direction
    end
    
    %% Constructor
    methods
        function this = LinOpDownsample(sz,df)
            this.name ='LinOpDownsample';
            assert(isequal(size(sz),size(df)),'Parameters sz and d must have the same size');
            assert(~any(mod(sz,df)),'Sizes in sz must be multiples of downsampling factors in df');  
            this.sizein=sz;
            this.df=df;
            this.sizeout=this.sizein./this.df;
            this.sel=cell(length(sz),1);
            for ii=1:length(sz)
                this.sel{ii}=1:this.df(ii):this.sizein(ii);
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
            y = zeros(this.sizein);
            y(this.sel{:}) = x;
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d,%d]',this.sizein);
            y = zeros(this.sizein);
            y(this.sel{:}) = x(this.sel{:});            
        end
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d, %d,%d,%d]',this.sizeout);
            y = x;
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            w=zeros(this.sizein);w(this.sel{:})=1;
            M=LinOpDiag([],w);
        end
    end
end
