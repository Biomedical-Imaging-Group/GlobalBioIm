classdef CostRectangle < CostIndicator
    %% CostRectangle :  Rectangle Indicator function
    %  Matlab Inverse Problems Library
    %
    % Obj = CostRectangle(xmin, xmax,H,y):
    % Implement the indicator function of the retangle set defined by xmin and xmax for the
    % $$ \phi(x) = 0 \textrm{ if } xmin <= (H x-y) <= xmax textrm{ and }  +\inf \textrm{ otherwise } $$
    % if x is complex
    % $$ \phi(x) = 0 \textrm{ if } (real(xmin) <= real(H x-y) <= real(xmax)
    %               && imag(xmin) <= imag(H x-y) <= imag(xmax))
    % textrm{ and }  +\inf \textrm{ otherwise } $$
    %
    %% Properties
    % * |xmin|         - minimum value (default -inf + 0i)
    % * |xmax|         - minimum value (default +inf + 0i)
    % * |iscomplex|    - true if acts on complex numbers 
    %
    %%
    %     Copyright (C) 2017 F. Soulez ferreol.soulez@epfl.ch
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
        xmin=-inf;  % lower bound (default -inf + 0i)
        xmax=+inf;  % higher bound  (default +inf + 0i)
        iscomplex = false;  % complex number flag
    end
    
    methods
        function this = CostRectangle(xmin, xmax,H,y)
            
            this.name='Cost  Rectangle';
            
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
                xmax=[];
            end
            if nargin==0
                xmin =-inf;
            end
            
            if isempty(xmax)
                xmax = +inf;
            end
            
            if((~isreal(xmin))||(~isreal(xmax)))
                this.iscomplex = true;
            end
            this.xmin = xmin;
            this.xmax = xmax;
            
        end
        
        
        function z = prox(this,x,~) % apply the operator
            
            if this.H.isinvertible
                
                if(isscalar(this.y)&&(this.y==0))
                    if this.iscomplex
                        res = this.H.apply(x);
                        tmp = min(max(real(res), real(this.xmin)),real(this.xmax)) + 1i.* min(max(imag(res), imag(this.xmin)),imag(this.xmax));
                        z = this.H.inverse(tmp);
                    else
                        z = this.H.inverse(min(max(real(this.H.apply(x)), this.xmin),this.xmax));
                    end
                else
                    if this.iscomplex
                        res = this.H.apply(x)-this.y;
                        tmp = min(max(real(res), real(this.xmin)),real(this.xmax)) + 1i.* min(max(imag(res), imag(this.xmin)),imag(this.xmax));
                        z = this.H.inverse(tmp+this.y);
                    else
                        z = this.H.inverse(min(max(real(this.H.apply(x)-this.y), this.xmin),this.xmax)+this.y);
                    end
                end
            else
                error('Prox not implemented');
            end
            
            
            
        end
        function cost = eval(this,x) % get the function cost
            cost = 0;
            
            if(isscalar(this.y)&&(this.y==0))
                res=this.H.apply(x);
            else
                res=this.H.apply(x)-this.y;
            end
            
            
            if this.iscomplex
                if(any(real(res(:))< real(this.xmin(:))))
                    cost= +inf;
                    return;
                end
                if(any(real(res(:))> real(this.xmax(:))))
                    cost= +inf;
                    return;
                end
                
                if(any(imag(res(:))> imag(this.xmax(:))))
                    cost= +inf;
                    return;
                end
                
                if(any(imag(res(:))< imag(this.xmin(:))))
                    cost= +inf;
                    return;
                end
                
            else
                if(~isreal(res))
                    if  any(imag(res(:)))
                        cost= +inf;
                    else
                        res = real(res);
                    end
                end
                               
                if(any(res(:)< this.xmin(:)))
                    cost= +inf;
                    return;
                end
                
                if(any(res(:)> this.xmax(:)))
                    cost= +inf;
                    return;
                end
                
            end
        end
    end
end
