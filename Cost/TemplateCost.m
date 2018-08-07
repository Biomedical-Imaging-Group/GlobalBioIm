classdef TemplateCost < Cost
	% TODO: Put here the description of your Cost
    % $$ C(\\mathrm{x}) = TODO $$
	%
	% :param parName: DESCRIPTION
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % **Note**: YOU CAN PUT A NOTE HERE
    %
    % **References**
    %
    % [1] Ref1 ...
    %
    % **Example** C=...
    %
    % See also :class:`Map` :class:`Cost`
    
    %%    Copyright (C) 
    %     TODO YEAR NAME EMAIL
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
 

    %% Properties
    % - Public
    properties (SetObservable, AbortSet)
        myVar
        % TODO : Set here new public properties
    end
    % - Readable
    properties (SetAccess = protected,GetAccess = public)
        % TODO : Set here new readable properties (read only)
    end
    % - Protected
    properties (SetAccess = protected,GetAccess = protected)
        % TODO : Set here new protected properties
    end
    
    %% Constructor
    methods     	
        function this = TemplateCost(sz,y)
        	% Default values
            if nargin<2, y=0; end
                % TODO : add others default values
            % Call superclass constructor
            this@Cost(sz,y);
            % Set properties
            this.name='TemplateCost';
            this.isConvex= ???;
            this.isDifferentiable=???;
            this.isSeparable=???;
                % TODO : set new defined properties 
            % Initialize
            this.initialize('TemplateCost'); % Call the method initialize with class name as attribute (du NOT put this.name BUT 'TemplateCost')

            % IMPORTANT : No computations in the constructor, only
            % affectations. Computations should be done in the method
            % updateProp (see below) in order to ensure a proper update 
            % when public properties are modified.
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`Cost`
            
            % Call superclass method
            updateProp@Cost(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'myVar') ||  strcmp(prop,'all')
                % TODO : code which has to be executed at each
                % modification of myVar
                % e.g. recomputation of the lipschitz constant
                this.lip= ???;
            end
            % TODO : other public properties
            
            % IMPORTANT : for each if statement ADD " || strcmp(prop,'all') "
            % in order to ensure that all computation are done at object
            % construction.
        end
    end
    
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            
            % TODO : IMPLEMENT IF APPLICABLE
        end
        function g=applyGrad_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            
            % TODO : IMPLEMENT IF APPLICABLE
        end
        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            
            % TODO : IMPLEMENT IF APPLICABLE
        end    
        
        % TODO : IMPLEMENT ANY USEFUL METHOD FROM INHERITED FROM MAP OR
        % COST CLASSES
    end
end
