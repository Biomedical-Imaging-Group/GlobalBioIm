classdef TestCvgStepRelative  < TestCvg
    % TestCvgStepRelative stops the optimization when the step  is below
    % the value STEPRELATIVETOL
    %
    % :param stepRelativeTol:  relative tolerance on step
    %
    % **Example** CvOpti=TestCvgStepRelative(stepRelativeTol )
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
    
    properties (SetAccess = public,GetAccess = public)
        stepRelativeTol=1e-5;      % stopping criteria tolerance on the relative difference btw two successive iterates
    end
    methods
        %% Constructor
        function this=TestCvgStepRelative( stepRelativeTol)
            this.name = 'TestCvgStepRelative';
            assert(isscalar(stepRelativeTol),'stepRelativeTol must be scalar');
            this.stepRelativeTol =stepRelativeTol;
            this.needxold = true;
        end
        %% Update method
        function stop = testConvergence(this,opti)
            % Tests algorithm convergence from the relative difference between two successive value of the step function
            %
            % :return: boolean true if
            %          $$ \\frac{\\| \\mathrm{x}^{k} - \\mathrm{x}^{k-1}\\|}{\\|\\mathrm{x}^{k-1}\\|} < \\text{stepRelativeTol}.$$
            
            stop = false;
            if ~isempty(opti.xold)
                r=opti.xopt-opti.xold;
                xdiff=norm(r(:))/(norm(opti.xold(:))+eps);
                if( xdiff < this.stepRelativeTol)
                    stop  =true;
                    endingMessage = [this.name,': Step variation between two successive iterates below the relative tolerance : ',num2str(xdiff),' < ',num2str(this.stepRelativeTol)];
                    opti.endingMessage = endingMessage;
                end
            end
        end
    end
end
