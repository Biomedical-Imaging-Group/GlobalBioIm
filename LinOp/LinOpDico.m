classdef LinOpDico <  LinOp
    % LinOpDico : Dictionnary Linear  Operator
    %
    % Defines a Dictionnary linear operator containing nbA atoms of size
    % szA according to [szA,nbA]=size(D) (i.e the last dimension of D is
    % used to index the different atoms).
    %
    % :param D: dictionary.
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** Dico=LinOpDico(D) 
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %%    Copyright (C) 
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

    properties (SetAccess = protected,GetAccess = public)
       D;         % Dictionnary
       ndms;      % number of dimensions of the dico
       numelAtom; % number of elements in one atom
    end
    
    %% Constructor
    methods
        function this = LinOpDico(D) 
            this.name ='LinOpDico';   
            sz=size(D);
            this.ndms=length(sz);
            this.sizein=[ones(1,this.ndms-1),sz(end)]; 
            this.sizeout = sz(1:end-1);   
            this.numelAtom=prod(this.sizeout);
            % Note: actually this operator can be invertible for some given
            % dictionnary D, but the inverse is not implemented 
            this.D=D;
            this.isDifferentiable=true;
            if this.numelAtom<1000 && sz(end)<1000  % the norm is computed if the dictionary is of reasonable size
                this.norm=norm(reshape(this.D,[this.numelAtom,sz(end)]));
            end
		end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        %todo: size problem here?
        function y = apply_(this,x)   
            % Reimplemented from parent class :class:`LinOp`.
            if cmpSize(size(x),this.sizein)
                y=sum(bsxfun(@times,this.D,x),3);
            else
                assert(isvector(x) && this.sizein(end)==length(x),'x must be a vector of size the number of dictionnary atoms');
                y=sum(bsxfun(@times,this.D,reshape(x,this.sizein)),3);
            end
        end		
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=sum(reshape(bsxfun(@times,this.D,x),[this.numelAtom,ones_(this.ndms-2),this.sizein(end)]),1);
            y=y(:);
        end
    end
end