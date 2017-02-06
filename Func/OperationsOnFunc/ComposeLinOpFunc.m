classdef ComposeLinOpFunc < Func
    %% ComposeLinOpFunc : Compose a Functional with a linear operator
    %  Matlab Inverse Problems Library
    %
    % -- Example
    % G =  ComposeLinOpFunc(F,Hcomp)
    % where F is a FUNC object and Hcomp a LINOP one
    %
    % Please refer to the FUNC superclass for general documentation about
    % functional class
    % See also Func
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
 
    % Protected Set and public Read properties     
    properties (SetAccess = protected,GetAccess = public)
		F;      % Functional
    end
    
    methods 
    	%% Constructor
        function this = ComposeLinOpFunc(F,Hcomp)
            this.name='ComposeLinOpFunc';
            this.F=F;
            if ~isa(F.H,'LinOpIdentity'), assert(isequal(Hcomp.sizeout,F.sizein),'sizeout of Hcomp must match with sizein of F'); end
			this.isconvex= F.isconvex; 
			this.set_H(Hcomp);
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			y=this.F.eval(this.H.Apply(x));
        end
    end
end
