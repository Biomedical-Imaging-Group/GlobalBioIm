classdef CostRing < Cost
    %% Complex Ring indicator 
    %  Matlab Inverse Problems Library
    %
    % Obj = CostRing(inner, outer,H,y):
    % Implement the proximity operator on the complex ring of inner radius INNER and outer radius OUTER
    % $$ \phi(x) = 0 \textrm{ if } inner<=||x||_2<=outer  textrm{ and }  +\inf \textrm{ otherwise } $$
    % Following Matlab angle function if x=0 then Circle(x) = 1
    %
    %% Properties
    %
    %%
    
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
    
    properties (SetAccess = protected,GetAccess = public)
        inner=0;
        outer = 1;
    end
    
    methods
        function this = CostRing(inner, outer, H,y)
            this.name='Cost Ring';
             % -- Set entries
            if nargin<4
                y=0;
            end
            if nargin<3
                H=[];
            end
            set_y(this,y);
            set_H(this,H);
            
            if nargin<2
                inner=0;
            end 
            if nargin==0
                outer =1;
            end
            
            
                if isscalar(inner)
                    this.inner = inner;
                else
                    if isnumeric(inner)
                        this.inner = inner; 
                        if ~isempty(this.H.sizeout)
                            assert(isequal(this.H.sizeout,size(inner)), 'inner must be equal to H.sizeout');
                        else
                            this.sizein = size(inner);
                        end
                    else
                        error('C should be numeric');
                    end
                end
                
                if isscalar(outer)
                    this.outer = outer;
                else
                    if isnumeric(outer)
                        this.outer = outer; 
                        if ~isempty(this.H.sizeout)
                            assert(isequal(this.H.sizeout,size(outer)), 'outer must be equal to H.sizeout');
                        else
                            this.sizein = size(outer);
                        end
                    else
                        error('C should be numeric');
                    end
                end
        end
        function z = prox(this,x,~) % apply the operator
            assert(isnumeric(x),'x must be numeric');
            
          if this.H.isinvertible
              tmp =this.H.apply(x) - this.y;
            z = this.H.inverse(max(min(this.outer,abs(tmp)),this.inner) .* exp(1i .* angle(tmp))+this.y);
            else
        		error('Prox not implemented');
          end
            
            
            
        end
    end
end
