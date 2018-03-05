classdef TestCvg  < handle
    % TestCvg class  monitor convergence criterion during
    % optimization
    %
    % :param verbose: if true will display a message before stopping the algorithm.
    %
    % At each iterations of an optimization algorithm (see :class:`Opti` generic class),
    % the update method of an :class:`TestCvg` object will be executed in order to acheive user 
    % defined computations, e.g.,
    %
    % **Example** CvOpti=TestCvg(verbose) 
    %
    % **Important** The update method should have an unique imput that is the :class:`Opti` object in order to 
    % be generic for all Optimization routines. Hence the update method has access (in reading mode) 
    % to all the properties of :class:`Opti` objects.
    %
    % See also :class:`Opti`

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
        name = 'TestCvg'% name of the optimization algorithm
        verbose = true;   
    end
    properties (SetAccess = public,GetAccess = public)
    end
    
    methods
    	%% Constructor
        function this=TestCvg(verbose) 
            this.verbose = verbose;
        end
        %% Update method
        function stop = testConvergence(this,opti)            % Tests algorithm convergence from the relative difference between two successive iterates
          stop = false;
            
        end
    end
end
