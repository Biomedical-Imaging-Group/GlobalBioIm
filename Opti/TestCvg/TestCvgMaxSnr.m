classdef TestCvgMaxSnr  < TestCvg
    % TestCvgMaxSnr stops the optimization when the step  is below
    % the value MaxSnrTOL
    %
    % :param MaxSnrTol:  relative tolerance on step
    %
    % **Example** CvOpti=TestCvgMaxSnr(MaxSnrTol )
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
        MaxSnrTol=1e-5;      % stopping criteria tolerance on the relative difference btw two successive iterates
        ref;                 % reference signal
    end
    properties (SetAccess = protected, GetAccess = protected)
        normRef      % norm of the reference signal
        sizeRef      % Size of the reference signal
    end
    methods
        %% Constructor
        function this=TestCvgMaxSnr( MaxSnrTol,ref)
            this.name = 'TestCvgMaxSnr';
            assert(isscalar(MaxSnrTol),'MaxSnrTol must be scalar');
            this.MaxSnrTol =MaxSnrTol;
            this.ref = ref;
            this.normRef=norm(ref(:));
        end
        %% Update method
        function stop = testConvergence(this,opti)
            % Tests algorithm convergence using the SNR with
            % :returns stop: boolean true if
            % $$ \\log\\left(\\frac{\\| \\mathrm{ref}\\|}{\\| \\mathrm{x}^{k} - \\mathrm{ref}\\|}\\right)< \\text{MaxSnrTol}.$$
            
            stop = false;
            assert(checkSize(this.sizeRef, size(opti.xopt)), ' Reference signal and x should be conformable');
            if ~isempty(this.oldStep)
                xsnr=20*log10(this.normRef/norm(this.ref(:)-opti.xopt(:)));
                if( xsnr < this.MaxSnrTol)
                    stop  =true;
                    endingMessage = [this.name,': SNR below the  tolerance : ',num2str(xsnr),' < ',num2str(this.MaxSnrTol)];
                    opti.endingMessage = endingMessage;
                end
            end
        end
    end
end
