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
    
    %% Properties
    % - Public 
    properties (SetObservable, AbortSet)
        mapsCell;    % Cell of summed Maps
        alpha;       % Correcponding coefficients
    end
    % - Readable
    properties (SetAccess = protected,GetAccess = public)
        numMaps;     % Number of Maps
    end
    
    %% Constructor
    methods
        function this = MapSummation(Maps,alpha)
            % Default values
            if nargin == 1, alpha = 1; end
            % Set properties
            this.name ='MapSummation';
            this.sizein = Maps{1}.sizein;
            this.sizeout = Maps{1}.sizeout;
            this.alpha = alpha;
            this.mapsCell = Maps;
            this.isInvertible=false;
            % Initialize
            this.initialize('MapSummation');
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`Map`
            
            % Call superclass method
            updateProp@Map(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'mapsCell') ||  strcmp(prop,'all')
                allMaps = all( cellfun(@(x)(isa(x, 'Map')), this.mapsCell) );
                assert(iscell(this.mapsCell) && allMaps, 'Property mapsCell should be a cell array Map');
                % If any of the inputs is itself a sum, expand it
                eMaps = []; % expanded LinOp list
                newAlphas = [];
                for i = 1:length(this.mapsCell)
                    if isa(this.mapsCell{i}, 'MapSummation') && ~isa(this.mapsCell{i},'CostPartialSummation')
                        eMaps = [eMaps this.mapsCell{i}.mapsCell];
                        newAlphas = [newAlphas  this.mapsCell{i}.alpha*this.alpha(i)];
                    else
                        eMaps = [eMaps this.mapsCell(i)];
                        newAlphas = [newAlphas this.alpha(i)];
                    end
                end
                this.alpha = newAlphas;
                this.mapsCell = eMaps;
                this.numMaps = length(this.mapsCell);
                this.isDifferentiable= this.mapsCell{1}.isDifferentiable;
                for n =2:this.numMaps
                    assert(cmpSize(this.sizein,this.mapsCell{n}.sizein),'In property mapsCell %d-th input does not have consistent  sizein', n) ;
                    assert(cmpSize(this.sizeout,this.mapsCell{n}.sizeout),'%In property mapsCell d-th input does not have the consistent sizeout ', n);
                    this.isDifferentiable= this.mapsCell{n}.isDifferentiable && this.isDifferentiable;
                end
            end
            if strcmp(prop,'alpha') ||  strcmp(prop,'all')
                assert(isnumeric(this.alpha) && (isscalar(this.alpha) || (isvector(this.alpha) && (isempty(this.numMaps) || numel(this.alpha)== this.numMaps))),...
                    'second input should be a scalar or an array of scalars of the same size as the first input');
                if  isscalar(this.alpha)
                    this.alpha = repmat(this.alpha, 1, this.numMaps) ;
                end
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
            x = zeros_(this.sizeout);
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
end

