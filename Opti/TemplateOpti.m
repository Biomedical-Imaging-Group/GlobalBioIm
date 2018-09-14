classdef TemplateOpti < Opti
    % OptiNAME which minimizes :class:`Cost`of the form
    % $$ C(\\mathrm{x}) = TODO $$
    % 
    % :param parNAME: DESCRIPTION
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
    % **Note**: YOU CAN PUT A NOTE HERE
    %
    % **References**
    %
    % [1] Ref1 ...
    %
    % **Example** opti=...
    %
    % See also :class:`Opti`, :class:`OutputOpti`, :class:`Cost`
    
    %%    Copyright (C) TODO YEAR
    %     TODO NAME EMAIL
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
        function this=TemplateOpti(~)
            % Default values
            if nargin<2, end
                % TODO : add default values if any
            % Set properties
            this.name='TemplateOpti';
            this.cost=????;
            this.OutOp=?????;
                % TODO : set new defined properties
            % Initialize
            this.initObject('TemplateOpti'); % Call the method initObject with class name as attribute (du NOT put this.name BUT 'TemplateOpti')
            
            % IMPORTANT : No computations in the constructor, only
            % affectations. Computations should be done in the method
            % updateProp (see below) in order to ensure a proper update
            % when public properties are modified.
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`Opti`
            
            % Call superclass method
            updateProp@Opti(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'myVar') ||  strcmp(prop,'all')
                % TODO : code which has to be executed at each
                % modification of myVar
            end
            % TODO : other public properties
            
            % IMPORTANT : for each if statement ADD " || strcmp(prop,'all') "
            % in order to ensure that all computation are done at object
            % construction.
        end
    end
    
    %% Methods for optimization
    methods
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
        end
        function doIteration(this)
            % Reimplementation from :class:`Opti`.
        end
        function updateParams(this)
            % Updates the parameters of the algorithm at each iteration
            % (default: no update). This method can be overloaded to makes
            % some parameters varying during iterations (e.g. descent step,
            % lagrangian parameters...)
        end
    end
end
