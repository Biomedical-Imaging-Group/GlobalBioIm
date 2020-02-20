classdef LinOpSummation < MapSummation &  LinOp
    % LinOpSummation: Sum of linear operators
    % $$ \\mathrm{H}(\\mathrm{x}) = \\sum_i \\alpha_i \\mathrm{H}_i(\\mathrm{x}) $$
    %
    % :param LinOps:  cell of :class:`LinOp`
    % :param alpha:  array of coefficients
    %
    % **Example** L=LinOpSummation(LinOps,alpha)
    %
    % See also :class:`Map`, :class:`LinOp`, :class:`MapOpSummation`
    
    %%    Copyright (C) 2017
    %     F. Soulez ferreol.soulez@epfl.ch
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
    
    %% Constructor
    methods
        function this = LinOpSummation(LinOps,alpha)
            this@MapSummation(LinOps,alpha); 
            this.name ='LinOpSummation';			
			allLinOps = all( cellfun(@(x)(isa(x, 'LinOp')), LinOps) );
			assert(iscell(LinOps) && allLinOps, 'First input should be a cell array LinOp'); 
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyAdjoint_(this,x)
    % - applyHtH_(this,x) 
    % - applyHHt_(this,y) 
    % - makeAdjoint_(this)
    methods (Access = protected)      		
        function x = applyAdjoint_(this,y) 
            % Reimplemented from :class:`LinOp` 
            x =  zeros_(this.sizein);
            for n = 1:this.numMaps
                x = x + this.alpha(n) .* this.mapsCell{n}(1).applyAdjoint(y);
            end
        end	
        function M = makeAdjoint_(this)
            % Reimplemented from :class:`LinOp`
            adjointCell = cellfun(@(x)x',this.mapsCell,'UniformOutput',false);
            M=LinOpSummation(adjointCell,this.alpha);
        end
        % the function reimplementations below is needed because of
        % the multiple inheritance to specifies which method to use from
        % parent classes
        function M = makeComposition_(this,G)
            M=makeComposition_@MapSummation(this,G);
        end
        function M = plus_(this,G)
            % Reimplemented from :class:`LinOp` 
            
            if isa(G,'LinOp')
                M=[];
                if isa(G,'LinOpSummation')
                    M=this+G.mapsCell{1};
                    for ii=2:G.numMaps
                        M=M+G.mapsCell{ii};
                    end
                else
                    % Find elements of the same type
                    ind=find(strcmp(G.name,cellfun(@(T) T.name,this.mapsCell,'UniformOutput',false)));
                    if length(ind)==1
                        M=G + this.mapsCell{ind};
                        if ~isa(M,'LinOpSummation')
                            % To avoid infinite loops (normally it should never goes in the else because the sum of
                            % two LinOp of the same type can always be simplified. If not the sum_ method of the corresponding
                            % LinOp has to be implemented properly).
                            for ii=1:this.numMaps
                                if ii~=ind
                                    M= M+this.mapsCell{ii};
                                end
                            end
                        else
                            M=[];
                        end
                    end
                end
                if isempty(M)
                    M=LinOpSummation({this,G},[1,1]);
                end
            else
                M = MapSummation({this,G},[1,1]);
            end
        end
    end
    
    methods (Access = protected)
        %% Copy
      function this = copyElement(obj)
          this = copyElement@MapSummation(obj);
          this.memoCache.applyAdjointInverse=struct('in', [], 'out', []);
          this.memoCache.applyAdjoint=struct('in', [], 'out', []);
          this.memoCache.applyHtH=struct('in', [], 'out', []);
          this.memoCache.applyHHt=struct('in', [], 'out', []);
      end
    end
end

