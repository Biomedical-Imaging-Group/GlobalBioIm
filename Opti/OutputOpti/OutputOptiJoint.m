classdef OutputOptiJoint < handle
    % OutputOptiJoint class for Joint algorithms displayings
    %
    % :param OutX: OutputOpti object associated to X
    % :param OutZ: OutputOpti object associated to Z
    %
    % **Note** This class only perform displaying. Saving of
    % iterates/cost/snr evolution are done inside OutX and OutZ.
    %
    % See also :class:`OptiJoint` :class:`OutputOpti`

    %%    Copyright (C) 2018 
    %     E. Soubies emmanuel.soubies@epfl.ch
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
    % - Protected
    properties (Access=protected)
        OutX  % OutputOpti for X
        OutZ  % OutputOpti for Z
    end
    
    %% Constructor, Initialization and Update method
    methods
        function this=OutputOptiJoint(OutX,OutZ)
            this.OutX=OutX;
            this.OutZ=OutZ;
        end
        function init(this)
            % Initialize the arrays and counters.
        	this.OutX.init();
            this.OutZ.init();
        end
        function update(this,opti)
            % Computes SNR, cost and display evolution.
            str=sprintf('Iter: (Out. %4i, Cum. %4i)',opti.niter,opti.cumIter);
            if this.OutX.computecost
                if this.OutX.iternum(end)==opti.cumIter  && ~isnan(this.OutX.evolcost(end))
                    str=sprintf('%s | Cost X: %+4.4e',str,this.OutX.evolcost(end));
                else
                    str=sprintf('%s | Cost X:    -----   ',str);
                end
            end
            if this.OutZ.computecost
                if this.OutZ.iternum(end)==opti.cumIter && ~isnan(this.OutZ.evolcost(end))
                    str=sprintf('%s | Cost Z: %+4.4e',str,this.OutZ.evolcost(end));
                else
                    str=sprintf('%s | Cost Z:    -----   ',str);
                end
            end
            if this.OutX.isgt
                if this.OutX.iternum(end)==opti.cumIter  && ~isnan(this.OutX.evolsnr(end))
                    str=sprintf('%s | SNR X: %+4.4e dB',str,this.OutX.evolsnr(end));
                else
                    str=sprintf('%s | SNR X:     -----     ',str);
                end
            end
            if this.OutZ.isgt
                if this.OutZ.iternum(end)==opti.cumIter  && ~isnan(this.OutZ.evolsnr(end))
                    str=sprintf('%s | SNR Z: %+4.4e dB',str,this.OutZ.evolsnr(end));
                else
                    str=sprintf('%s | SNR Z:     -----     ',str);
                end
            end
            if opti.verbose 
                disp(str);
            end
        end
    end

end