classdef CostL1 < Cost
    % L1 norm cost function
    % $$C(x) := \\|\\mathrm{Hx} - \\mathrm{y}\\|_1 $$
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param nonneg: boolean (varargin parameter) to combine a nonnegativity constraint to the cost (default false).
	%
    % See also :class:`Cost` :class:`LinOp`    

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
    end
    
    methods
        function this = CostL1(H,y,varargin)
            this.name='Cost L1';
            this.isconvex=true;
            % -- Set entries
            if nargin<2
                y=0;
            end
            if nargin<1
                H=[];
            end
            set_y(this,y);
            set_H(this,H);            
            
            for c=1:length(varargin)
                switch varargin{c}
                    case('NonNegativity')
                        this.nonneg = true;
                end
            end
            
        end
        
        function z = prox(this,x,alpha) 
        	% Reimplemented from parent class :class:`Cost` if the operator :attr:`H`  is invertible.
        	
            assert(isscalar(alpha),'alpha must be a scalar');           
            if this.H.isInvertible
                 tmp = this.H.apply(x)-this.y ;
                if this.nonneg
                    z =   this.H.inverse(max(tmp- alpha,0)+this.y);
                else
                    z =  this.H.inverse(sign(tmp) .* max( abs( tmp) - alpha,0)+this.y);
                end
            else
                error('Prox not implemented');
            end
        end
        function f=eval(this,x)
        	% Reimplemented from parent class :class:`Cost`.
        	
            if ~this.nonneg
                 tmp = this.H.apply(x)-this.y;
                 f=sum(abs(tmp(:)));
            else
                error('eval L1 not implemented');
            end                
        end       
    end
end
