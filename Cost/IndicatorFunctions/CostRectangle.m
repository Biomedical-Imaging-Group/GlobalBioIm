classdef CostRectangle < CostIndicator
    % CostRectangle:  Rectangle Indicator function
    % $$ C(x) = \\left\\lbrace \\begin{array}[l]
    % \\text{0~if } \\mathrm{real(xmin)} \\leq \\mathrm{imag(x-y)} \\leq \\mathrm{real(xmax)}
    % \\text{ and }  \\mathrm{imag(xmin)} \\leq \\mathrm{imag(x-y)} \\leq \\mathrm{imag(xmax)},  \\newline
    % + \\infty \\text{ otherwise.} \\end{array} \\right. $$
    %
    % :param xmin: minimum value (default -inf + 0i)
    % :param xmax: maximum value (default +inf + 0i)
    %
    % All attributes of parent class :class:`CostIndicator` are inherited 
    %
    % **Example** C=CostRectangle(sz,xmin,xmax,y)   
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostIndicator`

    %%    Copyright (C) 2017 
    %     F. Soulez ferreol.soulez@epfl.ch
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
        xmin=-inf;          % lower bound (default -inf + 0i)
        xmax=+inf;          % higher bound  (default +inf + 0i)         
    end
    properties (Access = protected)
        isComplex = false;  % complex number flag      
    end
    
    %% Constructor
    methods
        function this = CostRectangle(sz,xmin,xmax,y)   
            if nargin<4, y=0; end
            this@CostIndicator(sz,y);
            this.name='CostRectangle';                      
            if nargin<3, xmax=[]; end
            if nargin<2, xmin =-inf; end            
            if isempty(xmax), xmax = +inf; end           
            if((~isreal(xmin))||(~isreal(xmax)))
                this.isComplex = true;
            end
            this.xmin = xmin;
            this.xmax = xmax;
            this.isConvex=true;  
            this.isSeparable=true;
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,z,alpha)
    methods (Access = protected)
        function y = apply_(this,x) 
            % Reimplemented from parent class :class:`Cost`.  
            y = 0;          
            if(isscalar(this.y)&&(this.y==0))
                res=x;
            else
                res=x-this.y;
            end                  
            if this.isComplex
                if(any(real(res(:))< real(this.xmin(:))))
                    y= +inf;
                    return;
                end
                if(any(real(res(:))> real(this.xmax(:))))
                    y= +inf;
                    return;
                end               
                if(any(imag(res(:))> imag(this.xmax(:))))
                    y= +inf;
                    return;
                end               
                if(any(imag(res(:))< imag(this.xmin(:))))
                    y= +inf;
                    return;
                end               
            else
                if(~isreal(res))
                    if  any(imag(res(:)))
                        y= +inf;
                    else
                        res = real(res);
                    end
                end                              
                if(any(res(:)< this.xmin(:)))
                    y= +inf;
                    return;
                end               
                if(any(res(:)> this.xmax(:)))
                    y= +inf;
                    return;
                end            
            end
        end
        function y = applyProx_(this,x,~) 
            % Reimplemented from parent class :class:`Cost`.  
            if(isscalar(this.y)&&(this.y==0))
                if this.isComplex
                    y = min(max(real(x), real(this.xmin)),real(this.xmax)) + 1i.* min(max(imag(x), imag(this.xmin)),imag(this.xmax));
                else
                    y = min(max(real(x),this.xmin),this.xmax);
                end
            else
                if this.isComplex
                    res = x-this.y;
                    y = min(max(real(res), real(this.xmin)),real(this.xmax)) + 1i.* min(max(imag(res), imag(this.xmin)),imag(this.xmax));
                else
                    y = min(max(real(x-this.y), this.xmin),this.xmax);
                end
            end
        end       
    end
end
