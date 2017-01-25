classdef Func < handle
    %% Func : Functional generic class
    %  Matlab Inverse Problems Library
    %  The Func meta class implements generic methods for all functionals
    %
    % -- Properties
    % * |name|       - name of the function  
    % * |sizein|     - size of input space (kernel)
    % * |H|          - LinOp composed with the functional
    % * |lip|        - Lipschitz constant of the gradient (if known, otherwise -1)
    %
    % -- Methods
    % * |eval|       - evaluates the functional 
    % * |grad|       - evaluates the gradient of the functional 
    % * |prox|       - computes the proximity operator
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
    
    properties (SetAccess = protected,GetAccess = public)
        name = 'none'       % name of the functional
        H=LinOpIdentity();  % linear operator
        sizein;             % dimensions of the input vector space
        lip=-1;             % Lipschitz constant of the gradient
        % LinOp Infos
        isIdH=true;         % boolean (true if the linOp is identity)
    end
    
    methods
    	%% Evaluation of the Functional
        function eval(~,~) 
            error('Eval not implemented');
        end
        %% Gradient of the Functional
        function grad(~,~) 
            error('Prox not implemented');
        end
        %% Proximity operator of the functional
        function prox(~,~) 
            error('Prox not implemented');
        end
        %% Operator compose with a LinOp
        function v=	o(this,x)
        	error('operator compose (.o) not implemented');
        end
        %% Overload the operator +
        function y = plus(this,x)
            assert(isa(x,'Func'),'Addition of Func is only define with other Func');
            y = FuncSum({this,x});
		end
    end
end
