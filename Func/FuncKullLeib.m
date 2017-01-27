classdef FuncKullLeib < Func
    %% FuncKullLeib : Kullback–Leibler divergence
    %  Matlab Inverse Problems Library
    %
    % -- Description
    % Implements the Kullback–Leibler divergence
    % $$ \sum_n -d_n log((H*x)_n + bet) + (H*x)_n,
    % where H is a LinOp object (default LinOpIdentity), d are the data
    % and bet is a small scalar (default 1e-10) to smooth the function at zero.
    %
    % -- Example
    % F=FuncKullLeib(H,d,bet)
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
        data          % data 
		bet=1e-10;     % smoothing parameter
    end
    % Full protected properties 
    properties (SetAccess = protected,GetAccess = protected)
    	He1;     % Application of the adjoint of H to a vector of ones
    end
    
    methods 
    	%% Constructor
        function this = FuncKullLeib(data,H,bet)
            this.name='Func Kullback–Leibler';
            this.isconvex=true; 
            % -- Set entries
            if nargin==1 || isempty(H)
            	H=LinOpIdentity(); 
            end   
            if nargin==3, this.bet=bet; end
            this.data=data;
			this.set_H(H);
            assert( isequal(size(data),this.H.sizeout),'H sizeout and data size are not equal');
            tmp=ones(this.sizein);
            this.He1=this.H.Adjoint(tmp);
            % -- Compute Lipschitz constant of the gradient (if the norm of H is known)
            if this.H.norm>=0;
            	this.lip=max(this.data(:))/this.bet^2*this.H.norm^2;
            end
    	end
    	%% Evaluation of the Functional
        function y=eval(this,x)
			tmp=this.H.Apply(x);
			y=sum(-this.data(:).*log(tmp(:)+this.bet) + tmp(:));
        end
        %% Gradient of the Functional
        function g=grad(this,x)
			g= this.He1 - this.H.Adjoint(this.data./(this.H.Apply(x)+this.bet));
        end
        %% Proximity operator of the functional
        function y=prox(this,x,alpha)
        	y=[];
			if this.isIdH % if operator identity
				delta=(x-alpha-this.bet).^2+4*(x*this.bet + alpha*(this.data-this.bet));
				y=zeros(size(x));
				mask=delta>=0;
				y(mask)=0.5*(x(mask)-alpha-this.bet + sqrt(delta(mask)));
			end
			if isempty(y),error('Prox not implemented');end
        end
    end
end
