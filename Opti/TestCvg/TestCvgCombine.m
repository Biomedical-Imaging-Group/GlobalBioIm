classdef TestCvgCombine < TestCvg
    % TestCvgCombine:
    % Combine several tests 
    %
    %
    % At each iterations of an optimization algorithm (see :class:`Opti` generic class),
    % the update method of an :class:`TestCvg` object will be executed in order to acheive user
    % defined computations, e.g.,
    %
    % **Example** 
    % CvOpti=TestCvgCombine(A, B ,C) 
    % where A B and C are of TestCvg class
    %
    % CvTestVMLMB = TestCvgCombine('CostRelative',0.000001, 'CostAbsolute',10);
    % for simple test
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
        cvList = {};
        testNumber = 0;
    end
    properties (SetAccess = public,GetAccess = public)
    end
    
    methods
        %% Constructor
        function this=TestCvgCombine(varargin)
        this.name = 'TestCvgCombine';
            for n=1:nargin
                if isa(varargin{n},'TestCvg')
                    this.testNumber  = this.testNumber +1;
                    this.cvList{this.testNumber} = varargin{n};
                elseif ischar(varargin{n}) && exist(['TestCvg' varargin{n}],'class')
                    this.testNumber  = this.testNumber +1;
                    assert(n+1<=nargin,'the name of the TestCvg must be followed by a parameter');
                    this.cvList{this.testNumber} = eval(['TestCvg' varargin{n} '(' num2str(varargin{n+1}) ')']) ;
                    n = n+1;
                elseif isa(varargin{n},'TestCvgCombine')
                    this.cvList = {this.cvList, varargin{n}.cvList };
                    this.testNumber = this.testNumber + varargin{n}.testNumber;
                end
            end
        end
        %% Update method
        function stop = testConvergence(this,opti)            % Tests algorithm convergence from the relative difference between two successive iterates
            stop = false;
            for n=1:this.testNumber
                stop = this.cvList{n}.testConvergence(opti);
                if stop, break; end
            end
        end
    end
end

