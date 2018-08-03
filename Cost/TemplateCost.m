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
    properties (SetAccess = protected,GetAccess = protected)
        % TODO : Set here new readable properties (read only)
    end
    % - Private
    properties (SetAccess = protected,GetAccess = protected)
        % TODO : Set here new private properties
    end
    
    %% Constructor and handlePropEvents method
    methods     	
        function this = TemplateCost(sz,y)
        	% Default values
            if nargin<2, y=0; end
                % TODO : add others default values
            % Call superclass constructor
            this@Cost(sz,y);
            % Listeners to PostSet events
            addlistener(this,'myVar','PostSet',@this.handlePropEvents);
                % TODO : add a 'PostSet' listener for all public properties
            % Set properties
            this.name='CostNAME';
            this.isConvex= ???;
            this.isDifferentiable=???;
            this.isSeparable=???;
                % TODO : set new defined properties 
            % Listeners to modified events (for properties which are classes) 
            addlistener(this.myVar,'modified',@this.handleModifiedmyVar);
                % TODO : add a 'modified' listener for all public
                % properties which are classes of the library (LinOp, Cost,
                % Opti, Map)

            % IMPORTANT : No computations in the constructor, only
            % affectations. Computations should be done in the method
            % handlePropEvents (see below) in order to ensure a proper update 
            % when public properties are modified.
        end
        function handleModifiedmyVar(this,~,~) % Necessary for properties which are objects of the Library
            sourc.Name='myVar'; handlePropEvents(this,sourc);
        end
        % TODO : add one method like that for each 'modified' listener
        function handlePropEvents(this,src,~)
            % Reimplemented from superclasses :class:`Map` and :class:`MapComposition`
            switch src.Name
                case 'myVar'
                    % TODO : code which has to be executed at each
                    % modification of myVar
                    % e.g. recomputation of the lipschitz constant
                    this.lip= ???; 
                % TODO : other public properties
            end
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
