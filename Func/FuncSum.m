classdef FuncSum < Func
    %% FuncSum : Sum of Functionnals
    %  Matlab Inverse Problems Library
    %
    % -- Example
    % F = FuncSum(AFunc,alpha)
    % Sum the all Func contained in vector AFUNC weighted by ALPHA (default 1)
    % F  sum_n alpha(n) * AFunc(n)
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
        funcs      % cell of Func
        numfuncs   % number of Func
        alpha      % scalar factor (array)
    end
    
    methods 
    	%% Constructor
        function this = FuncSum(funcs,alpha)
            if nargin == 1, alpha = ones(size(funcs));end
            this.name='Func Sum';
            this.numfuncs = numel(funcs);
            assert(isnumeric(alpha)&& ( isscalar(alpha) || ( isvector(alpha) && (numel(alpha)== this.numfuncs))),'Second input should be a scalar or an array of scalar of the same size as the first input');
            allfuncs = all( cellfun(@(x)(isa(x, 'Func')), funcs) );
			assert(iscell(funcs) && allfuncs, 'First input should be a cell array Func');
			if  isscalar(alpha)
				this.alpha = repmat(alpha, 1, this.numfuncs) ;
			else
				this.alpha = alpha;
			end
			this.funcs = funcs;
			this.sizein =  this.funcs{1}.sizein;
			this.isconvex=funcs{1}.isconvex; 
            for n =2:this.numfuncs
            	if isempty(this.sizein), this.sizein=this.funcs{n}.sizein; end
                assert(isempty(this.funcs{n}.sizein)  || isequal(this.sizein,this.funcs{n}.sizein),'%d-th input does not have the right hand side size ', n) ;
                this.isconvex = this.isconvex & funcs{1}.isconvex; 
            end
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			y=this.alpha(1)*this.funcs{1}.eval(x);
			for n=2:this.numfuncs
				y=y+this.alpha(n)*this.funcs{n}.eval(x);
			end
        end
        %% Gradient of the Functional
        function g=grad(this,x)
			g=this.alpha(1)*this.funcs{1}.grad(x);
			for n=2:this.numfuncs
				g=g+this.alpha(n)*this.funcs{n}.grad(x);
			end
        end
    end
end
