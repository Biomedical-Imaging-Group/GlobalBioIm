classdef CostNonNeg < CostReals
    %% CostNonNeg : Non negativity indicator 
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implement the indicator over positive vector function:
    % $$ \phi(Hx) = 0 \textrm{ if } Hx \ge 0 textrm{ and }  +\inf \textrm{ otherwise } $$
    %
    % -- Example
    % F = CostNonNeg();
    %
    % Please refer to the COST superclass for general documentation about
    % functional class
    % See also Cost, CostIndicator, LinOp
	%
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
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
    
    methods 
    	%% Constructor
        function this = CostNonNeg(H,y) 
            if nargin<2
                y=0;
            end
            if nargin<1
                H=[];
            end
            set_y(this,y);
            set_H(this,H);

            this.name='Cost NonNegativity';
			 
            this.xmin=0;
			 this.xmax =+inf;
             this.iscomplex=false;
        end
    end
end
