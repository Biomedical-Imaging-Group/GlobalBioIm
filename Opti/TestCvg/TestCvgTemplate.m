classdef TestCvgTemplate < TestCvg
    % TEMPLATE for the creation of a TestCvg to monitor convergence
    %
    %
    % At each iterations of an optimization algorithm (see :class:`Opti` generic class),
    % the update method of an :class:`TestCvg` object will be executed in order to acheive user
    % defined computations, e.g.,
    %
    % **Example** CvOpti=TestCvgTemplate()
    %
    % **Important** The update method should have an unique imput that is the :class:`Opti` object in order to
    % be generic for all Optimization routines. Hence the update method has access (in reading mode)
    % to all the properties of :class:`Opti` objects.
    %
    % See also :class:`TestCvg`
    
    %%    Copyright (C) 2018
    %     F. Soulez ferreol.soulez@univ-lyon1.fr
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
    end
    properties (SetAccess = public,GetAccess = public)
    end
    
    methods
        %% Constructor
        function this=TestCvgTemplate()
            this.name = 'TestCvgTemplate';
        end
        %% Update method
        function stop = testConvergence(this,opti)      
            stop = false;
            
            if condition
                stop = true;
                endingMessage = [this.name,': reason to stop : '];
                opti.endingMessage = endingMessage;
            end
        end
    end
end
