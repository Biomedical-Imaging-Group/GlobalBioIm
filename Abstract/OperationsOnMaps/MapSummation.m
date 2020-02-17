classdef MapSummation < Map
    % MapSummation: Sum of Maps
    % $$ \\mathrm{H}(\\mathrm{x}) = \\sum_i \\alpha_i \\mathrm{H}_i(\\mathrm{x}) $$
    %
    % :param Maps:  cell of :class:`Map`
    % :param alpha:  array of coefficients
    %
    % **Example** H=MapSummation(Maps,alpha)
    %
    % See also :class:`Map`, :class:`LinOpSummation`
    
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
    
    properties(SetAccess = protected,GetAccess = public)
        mapsCell;    % Cell of summed Maps
        alpha;       % Correcponding coefficients
        numMaps;     % Number of Maps
    end
    %% Constructor
    methods
        function this = MapSummation(Maps,alpha)
            this.name ='MapSummation';
            if nargin == 1
				alpha = 1;
            end			
			this.numMaps = numel(Maps);
            % Check coefs alpha
			assert(isnumeric(alpha)&& ( isscalar(alpha) || ( isvector(alpha) && (numel(alpha)== this.numMaps))),...
                'second input should be a scalar or an array of scalars of the same size as the first input');
			if  isscalar(alpha)
				this.alpha = repmat(alpha, 1, this.numMaps) ;
			else
				this.alpha = alpha;
			end
			% Check type Maps
			allMaps = all( cellfun(@(x)(isa(x, 'Map')), Maps) );
			assert(iscell(Maps) && allMaps, 'First input should be a cell array Map');		
			% If any of the inputs is itself a sum, expand it
			eMaps = []; % expanded LinOp list
			newAlphas = [];
			for i = 1:length(Maps)
				if isa(Maps{i}, 'MapSummation') && ~isa(Maps{i},'CostPartialSummation')
					eMaps = [eMaps Maps{i}.mapsCell];
					newAlphas = [newAlphas  Maps{i}.alpha*this.alpha(i)];
				else
					eMaps = [eMaps Maps(i)];
					newAlphas = [newAlphas this.alpha(i)];
				end
			end
			this.alpha = newAlphas;
			this.mapsCell = eMaps;
			this.numMaps = length(this.mapsCell);
			% Set some properties
            this.isDifferentiable= this.mapsCell{1}.isDifferentiable;
            this.isInvertible=false;
            this.sizein = this.mapsCell{1}.sizein;
            this.sizeout = this.mapsCell{1}.sizeout;
            this.norm=0;
            for n=1:this.numMaps
                if this.mapsCell{n}.norm~=-1
                    this.norm=this.norm+abs(this.alpha(n))*this.mapsCell{n}.norm;
                else
                    this.norm=-1;
                    break
                end
            end
            for n =2:this.numMaps
                assert(cmpSize(this.sizein,this.mapsCell{n}.sizein),'%d-th input does not have consistent  sizein', n) ;
                assert(cmpSize(this.sizeout,this.mapsCell{n}.sizeout),'%d-th input does not have the consistent sizeout ', n);
                this.isDifferentiable= this.mapsCell{n}.isDifferentiable && this.isDifferentiable;
            end      
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - makeComposition_(this,G)
    methods (Access = protected)
        function y = apply_(this,x) 
            % Reimplemented from :class:`Map`   
            y = zeros_(this.sizeout);
            for n = 1:this.numMaps
                y = y + this.alpha(n) .* this.mapsCell{n}.apply(x);
            end
        end  
        function x = applyJacobianT_(this, y, v)
            % Reimplemented from :class:`Map`   
            x = zeros_(this.sizein);
            for n = 1:this.numMaps
                x = x + this.alpha(n) .* this.mapsCell{n}.applyJacobianT(y,v);
            end
        end     
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`Map`  
            M=this.alpha(1)*this.mapsCell{1}*G;
            for i=2:this.numMaps
                M=M+ this.alpha(i)*this.mapsCell{i}*G;
            end
        end
    end
    
    methods (Access = protected)
        %% Copy
      function this = copyElement(obj)
          this = copyElement@Map(obj);
            for n = 1:this.numMaps
                 this.mapsCell{n} = copy(obj.mapsCell{n});
            end
      end
    end
end

