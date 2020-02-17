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

    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		% TODO SET HERE NEW PROTECTED SET AND PUBLIC READ PROPERTIES IF NEEDED.
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
		% TODO SET HERE NEW FULLY PROTECTED PROPERTIES 
		% (E.G. INTERNAL VARIABLE USED TO AVOID MULTIPLE COMPUTATION)
    end
    
    %% Constructor
    methods
        function this = TemplateLinOp() %change the name of the constructor
            this.name ='';             
            this.isInvertible = ??; 
            this.sizein = ??;         
            this.sizeout = ??;         
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
        
        %% Deep copy: IMPLEMENT IF APPLICABLE
        function this = copyElement(obj)
            this = copyElement@LinOp(obj);
            % properties that are handle objects such as Map have to be copied explicitly to ensure
            % e.g. this.OBJECT = copy(obj.OBJECT)
        end
    end
end
