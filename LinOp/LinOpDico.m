classdef LinOpDico <  LinOp
    %% LinOpDico : Dictionnary Linear  Operator
    %  Matlab Linear Operator Library 
    %
    % -- Example
    % Op = LinOpDico(D)
    % Defines a Dictionnary linear operator containing nbA atoms of size
    % szA according to [szA,nbA]=size(D) (i.e the last dimension of D is
    % used to index the different atoms).
    %
    % Please refer to the LinOp superclass for documentation
    % See also LinOp   
    %
    %     Copyright (C) E. Soubies emmanuel.soubies@epfl.ch
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
       D;         % Dictionnary
       ndms;      % number of dimensions of the dico
       numelAtom; % number of elements in one atom
    end
    methods
        function this = LinOpDico(D) 
            this.name ='LinOp Dictionnary';   
            sz=size(D);
            this.ndms=length(sz);
            this.sizein=[ones(1,this.ndms-1),sz(end)]; 
            this.sizeout = sz(1:end-1);   
            this.numelAtom=prod(this.sizeout);
            % By default iscomplex and isinvertible are false
            if ~isreal(D)
                this.iscomplex= true;     
            end
            % Note: actually this operator can be invertible for some given
            % dictionnary D, but the inverse is not implemented 
            this.D=D;
            if this.numelAtom<1000 && sz(end)<1000  % the norm is computed if the dictionary is of reasonable size
                this.norm=norm(reshape(this.D,[this.numelAtom,sz(end)]));
            end
		end
	end
	
	methods (Access = protected)
        function y = apply_(this,x)   %todo: size problem here?
            if isequal(size(x),this.sizein)
                y=sum(bsxfun(@times,this.D,x),3);
            else
                assert(isvector(x) && this.sizein(end)==length(x),'x must be a vector of size the number of dictionnary atoms');
                y=sum(bsxfun(@times,this.D,reshape(x,this.sizein)),3);
            end
		end
		
        function y = adjoint_(this,x)
            y=sum(reshape(bsxfun(@times,this.D,x),[this.numelAtom,ones(this.ndms-2),this.sizein(end)]),1);
            y=y(:);
        end
    end
end