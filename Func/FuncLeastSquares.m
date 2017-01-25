classdef FuncLeastSquares < Func
    %% FuncLeastSquares : Least Squares functional
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implement the cost function for the weighted L2 norm
    % $$ \phi(x) = 1/2||Hx - d||^2_W $$
    %
    % -- Example
    % Obj = FuncLeastSquares(H,d,W);
    % where H is a LinOp object (default LinOpIdentity), d are the data and W
    % is a weight matrix (LinOp, default LinOpIdentity).
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
    
    properties (SetAccess = protected,GetAccess = public)
        W=[];       % weight matrix
        WplusWt=[]; % sum of W plus its transpose
        isW=false;  % boolean 
        data        % data 
        Hd          % application of the adjoint to data
        % LinOp Infos
        isConvH=false; % boolean (true if the linOp is convolution)
    end
    
    methods 
    	%% Constructor
        function this = FuncLeastSquares(data,H,wght)   
            if nargin==1
            	H=[];wght=[];
            	this.Hd=this.H.Adjoint(data);
            	this.sizein=this.H.sizein;
            else
            	if ~isempty(H), this.H=H; end       	
            	if strcmp(this.H.name,'LinOp Convolution'), this.isConvH=true; end;
            	if nargin==3
            		this.W=wght;
            		this.WplusWt=wght+wght';
            		this.isW=true;
            		this.Hd=this.H.Adjoint(this.WplusWt.Apply(data)); 
            	else	
            		this.Hd=this.H.Adjoint(data);
            	end
            	assert( isequal(size(data),this.H.sizeout),'H sizeout and data size are not equal');
            	this.sizein=this.H.sizein;
            end
            this.data=data;
            this.name='Func Least Squares';
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
        	r=this.H.Apply(x)-this.data;
        	if this.isW
        		wr=this.W.Apply(r);
        		y=0.5*dot(r(:),wr(:));
        	else
				y=0.5*norm(r(:))^2;
			end
        end
        %% Gradient of the Functional
        function g=grad(this,x)
        	if this.isW
				g = 0.5*(this.H.HtWH(x,this.WplusWt)-this.Hd);
			else
				g = this.H.HtH(x) - this.Hd;
			end
        end
        %% Operator compose with a LinOp
        function v=	o(this,x)
        	assert(isa(x,'LinOp'),' Composition of Func(.o) is only define with a LinOp');
        	v=this;
        	v.H=this.H*x;
        	if this.isW
        		v.Hd=v.H.Adjoint(v.WplusWt.Apply(data)); 
        	else
        		v.Hd=v.H.Adjoint(this.data);
        	end
        end
        %% Proximity operator of the functional
        function y=prox(this,x)

        end
    end
end
