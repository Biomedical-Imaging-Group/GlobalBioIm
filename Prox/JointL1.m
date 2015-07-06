classdef JointL1 < Prox
    %% JointL1 :  Proximity Operator for the JointL1 norm
    %  Matlab Inverse Problems Library
    %
    %
    % Obj = JointL1(index):
    % Implement the proximity operator for the jointL1 norm aka L21 aka
    % group LASSO aka structured L1 aka group L1 ...
    % $$ \phi(x) = \sum_k \sqrt( \sum_l x^2_{k,l} ) $$
    % INDEX with indicate on which dimensions will the inner sum (l is the
    % example)
    % The option 'NonNegativity' add the non negativity constraint (default
    % false)
    
    %
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
    %%
    
    properties (SetAccess = protected,GetAccess = public)
        nonneg = false
        index
    end
    
    methods
        function this = JointL1(index,varargin)
            this.name='JointL1';
            
            assert(isvector(index),'The index should be a conformable  to sz');
            this.index = index;
             
            for c=1:length(varargin)
                switch varargin{c}
                    case('NonNegativity')
                        this.nonneg = true;
                end
            end
            
        end
        function y = Apply(this,x,alpha) % Apply the operator
            assert(isscalar(alpha),'alpha must be a scalar');
            this.alpha = alpha;
            
            
            sz = size(x);
            ndms = length(sz);
             T = true(ndms,1);
             T(this.index)=false;
           kerdims = sz;
            kerdims(T)=1;
            imdims = sz;
            imdims(~T)=1;
            
            
            
             if this.nonneg
                x =  max(x,0);
            end
            index = this.index;
            sx = abs(x).^2;
            while ~isempty(index)
                sx = sum(sx,index(1));
                index = index(2:end);
            end
            sx = sqrt(sx);
            t = sx > this.alpha;
            b = zeros(size(sx));
            
            b(t) = 1-this.alpha./sx(t);
            y = reshape(repmat(reshape(b ,imdims),kerdims),sz).*x;
            
            
             
        end
    end
end
