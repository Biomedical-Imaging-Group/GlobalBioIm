classdef CostMixNorm12 < Cost
    %% CostMixNorm12 : Mixed norm 1-2
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implements the mixed norm ||Hx||_{1,2}:
    % $$ \sum_k \sqrt( \sum_l (Hx)_{k,l}^2 ) $$
    % where H is a LinOp object (default LinOpIdentity)
    % It corresponds to an l2-norm along some dimensions (idexed by l in the example)
    % combined by an l1-norm along the remaining dimensions (sum over k in the example)
    %
    % -- Example
    % F = CostMixNorm21(index,H);
    % INDEX with indicate on which dimensions will the inner sum (l is the example)
    %
    % -- Properties
    % * |index|     dimensions along which the l2-norm will be applied
    %
    % Please refer to the FUNC superclass for general documentation about
    % functional class
    % See also Cost
	%
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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
 
    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		index;    % dimensions along which the l2-norm will be applied
    end
    
    methods 
    	%% Constructor
        function this = CostMixNorm12(index,H)
            this.name='Cost MixNorm1-2';
			this.isconvex=true; 
			if nargin==1 || isempty(H)
            	H=LinOpIdentity(); 
            end     
            assert(isnumeric(index)&&isvector(index),'The index should be a vector of integers');
            this.index=index;
            this.set_H(H);
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			u=abs(this.H.Apply(x)).^2;
			% Computes the l2-norm along the dimensions given by index
			for n=1:length(this.index)
				u = sum(u,this.index(n));
			end
			u = sqrt(u);
			y=sum(u(:));
        end
        %% Proximity operator of the functional
        function y=prox(this,x,alpha)
        	assert(isscalar(alpha),'alpha must be a scalar');
        	y=[];
        	if isa(this.H,'LinOpIdentity')
            	sz = size(x);
				ndms = length(sz);
				T = true(ndms,1);
				T(this.index)=false;
				kerdims = sz; kerdims(T)=1;
				imdims = sz; imdims(~T)=1;
			
				% Computes the l2-norm along the dimensions given by index
				sx = abs(x).^2;
				for n=1:length(this.index)
					sx = sum(sx,this.index(n));
				end
				sx = sqrt(sx);
				
				% Computes the prox
				t = sx > alpha;
				b = zeros(size(sx));
			
				b(t) = 1-alpha./sx(t);
				y = reshape(repmat(reshape(b ,imdims),kerdims),sz).*x;
			end
			if isempty(y),error('Prox not implemented');end
			% result:
			% x(||x|| <= alpha) = 0
			% x(||x|| > alpha) = x(||x|| > alpha) - ...
			%      x(||x|| > alpha) / ||x|| * alpha
        end
    end
end
