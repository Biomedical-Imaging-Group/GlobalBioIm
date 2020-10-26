classdef OpEWfunc < Map
    % Element-wise function operator: compute the point-wise function of
    % the input
    % 
    % :param sz: input size
    % :param func: the function to apply element-wise
    % :param func_grad: the derivative of the function to apply element-wise (default: non differentiable)
    % :param func_inv: the inverse of the function to apply element-wise (default: non invertible)
    %
    % All attributes of parent class :class:`Map` are inherited. 
    %
    % **Example** Mod = OpEWfunc(sz, func, func_grad, func_inv)
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2018 
    %     Created: 09/28/2020 (mm/dd/yyyy)
    %     Anthony Berdeu (National Astronomical Research Institute
    %     of Thailand)
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
        func ;          % Function
        func_grad ;     % Function gradient
        func_inv ;      % Function inverse
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
    end
    
    %% Constructor
    methods
        function this = OpEWfunc(sz, func, func_grad, func_inv)
            this.name ='OpEWfunc';
            this.sizein=sz;
            this.sizeout=sz;
            
            this.func = func ;
            
            if nargin < 3 || isempty(func_grad)
                this.isDifferentiable=false ;
                this.func_grad = [] ;
            else
                this.isDifferentiable=true ;
                this.func_grad = func_grad ;
            end
            if nargin < 4 || isempty(func_inv)
                this.isInvertible=false ;
                this.func_inv = [] ;
            else
                this.isInvertible=true ;
                this.func_inv = func_inv ;
            end
        end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`Map`.
            y=this.func(x) ;
        end	
        function x = applyJacobianT_(this,y,v)
            % Reimplemented from parent class :class:`Map`.
            x=this.func_grad(v).*y ;
        end	
        function x = applyInverse_(this,y)
            % Reimplemented from parent class :class:`Map`.
            x=this.func_inv(y) ;
        end	
    end
end
