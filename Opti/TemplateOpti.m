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

    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		% TODO SET HERE NEW PROTECTED SET AND PUBLIC READ PROPERTIES
		% IF NEEDED.
		% EXAMPLE THE MINIMIZED FUNC
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
		% TODO SET HERE NEW FULLY PROTECTED PROPERTIES 
		% (E.G. INTERNAL VARIABLE USED TO AVOID MULTIPLE COMPUTATION)
    end
    % Full public properties
    properties
		% TODO SET FULLY PUBLIC PROPERTIES
    end
    
    methods
    	%% Constructor
    	function this=TemplateOpti(~)
    		% TODO SET THE INHERITED PROPERTIES
    		this.name='TemplateOpti';
    		this.cost=????;
    		this.OutOp=?????;
    		% TODO SET NEW DEFINED PROPERTIES
        end
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
