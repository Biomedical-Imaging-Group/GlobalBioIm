classdef L2 < Prox
    %% L2 : L2 norm Proximity Operator 
    %  Matlab Inverse Problems Library
    % 
    % Obj = L2(a,wght):
    % Implement the proximity operator for the weighted L2 norm
    % $$ \phi(x) = 1/2||x - a||^2_wght $$
    % 
   
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
        a     % 
        wght  % weight
    end
    
    methods 
        function this = L2(a,wght)
            this.name='L2';
            switch nargin
                case 0
                    a =0;
                    wght=1;
                case 1
                    wght=1;
            end
            if isscalar(a)
                this.a = a;
            else
                this.a = a;
                this.sz = size(a);
            end
            if isscalar(wght)
                this.wght = wght;
            else
                assert( (~isscalar(a))&&(isequal(size(wght),this.sz)),  'x does not have the right size: [%d, %d]',this.sz);
                this.wght = wght;
                this.sz = size(wght);
            end
        end
        function y = Apply(this,x,alpha) % Apply the operator
            assert(isscalar(alpha),'alpha must be a scalar');
            this.alpha = alpha;
            if ~isscalar(this.a)
               assert( isequal(size(x),this.sz),  'x does not have the right size: [%d, %d, %d,%d]',this.sz); 
            end
            y = 1./(this.alpha + this.wght) .* (this.alpha.* x + this.wght .* this.a);
        end
        function cost = FCost(this,x) 
           cost =   this.wght .*(abs(x - this.a)).^2;
           cost = 0.5 * sum(cost(:));
        end
    end
end
