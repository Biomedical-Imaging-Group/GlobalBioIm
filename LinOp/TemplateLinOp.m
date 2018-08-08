classdef TemplateLinOp <  LinOp
    % TODO: Put here the description of your LinOp
    % 
    % :param parName: DESCRIPTION
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note**: YOU CAN PUT A NOTE HERE
    %
    % **References**
    %
    % [1] Ref1 ...
    %
    % **Example** C=...
    %
    % See also :class:`Map` :class:`LinOp`
  
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
        function this = TemplateLinOp(sz,par1) %change the name of the constructor
            % Default values
            if nargin<2, par1=0; end
                % TODO : add others default values
            % Set properties
            this.name ='TemplateLinOp';             
            this.isInvertible = ??; 
            this.isDifferentiable=???;
            this.sizein = ??;         
            this.sizeout = ??;  
                % TODO : set new defined properties 
            % Initialize
            this.initialize('TemplateLinOp'); % Call the method initialize with class name as attribute (du NOT put this.name BUT 'TemplateLinOp')

            % IMPORTANT : No computations in the constructor, only
            % affectations. Computations should be done in the method
            % updateProp (see below) in order to ensure a proper update 
            % when public properties are modified.
		end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`LinOp`
            
            % Call superclass method
            updateProp@LinOp(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'myVar') ||  strcmp(prop,'all')
                % TODO : code which has to be executed at each
                % modification of myVar
                % e.g. recomputation of the lipschitz constant
                this.norm= ???;
            end
            % TODO : other public properties
            
            % IMPORTANT : for each if statement ADD " || strcmp(prop,'all') "
            % in order to ensure that all computation are done at object
            % construction.
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
	methods (Access = protected)
        function y = apply_(this,x)   
            % Reimplemented from parent class :class:`LinOp`.
            
            % TODO : IMPLEMENT IF APPLICABLE
        end		
        function y = applyAdjoint_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            % TODO : IMPLEMENT IF APPLICABLE
		end
        function y = applyHtH_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            % TODO : IMPLEMENT IF APPLICABLE
		end
        function y = applyHHt_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            % TODO : IMPLEMENT IF APPLICABLE
		end
        function y = applyInverse_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            % TODO : IMPLEMENT IF APPLICABLE
		end
		
        % TODO : IMPLEMENT ANY USEFUL METHOD FROM INHERITED FROM MAP OR
        % LINOP CLASSES
    end
end
