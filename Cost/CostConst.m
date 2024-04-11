classdef CostConst < Cost
    % CostConst: Constant cost whatever the input. Can be used to generate
    % a null cost in particular situations.
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param sizein: size of the input vector
    %
    % :param const: the output constant value (default 0)
    % (default true)
    %
    % **Example** C=CostConst(sz, const)
    %
    % **Example** C=CostMod2(sz)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`LinOp`
    
    %%    Copyright (C) 2018
    %     Created: 09/21/2018 (mm/dd/yyyy)
    %     Anthony Berdeu (Laboratoire Hubert Curien)
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
        const;      % the constant output value
    end
    
    %% Constructor
    methods     	
        function this = CostConst(sizein, const)
            
            %% Declaration
            this.name='CostConst';
            this.sizein = sizein ;
            
            this.isConvex=true;
            this.isDifferentiable=true;
            
            %% Checking inputs
            % const
            if nargin<2 || isempty(const)
                const = 0 ;
            end
            this.const = const ;
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        function y=apply_(this, ~)
            % Reimplemented from parent class :class:`Cost`.
            
            y = this.const ;
        end
        function g=applyGrad_(this,~)
            % Reimplemented from parent class :class:`Cost`.
            
            g = zeros(this.sizein) ;
        end
        function y=applyProx_(~,x,~)
          % Reimplemented from parent class :class:`Cost`.
             y = x;
        end
    end
end
