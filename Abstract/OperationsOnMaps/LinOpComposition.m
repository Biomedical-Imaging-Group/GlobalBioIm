classdef LinOpComposition < MapComposition & LinOp
    % LinOpComposition : Composition of LinOps
    % $$ \\mathrm{H}(\\mathrm{x}) = \\mathrm{H}_1 \\mathrm{H}_2\\mathrm{x} $$
    %
    % :param H1:  left hand side :class:`LinOp` (or a scalar)
    % :param H2:  right hand side :class:`LinOp`
    %
    % See also :class:`Map`, :class:`LinOp`, :class:`MapComposition`
    
    %%    Copyright (C) 2015
    %     F. Soulez ferreol.soulez@epfl.ch
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
		isHTH = 0;
		isHHt = 0;
    end
    
    %% Constructor
    methods
        function this = LinOpComposition(H1,H2)
            this@MapComposition(H1,H2);               
            assert((isa(H1,'LinOp') || isscalar(H1)),'H1 have to be a LinOp object or a scalar');
            assert(isa(H2,'LinOp'),'H2 have to be a LinOp');
            if (isnumeric(H1) && isscalar(H1)) && isa(H2,'LinOpComposition') && (isnumeric(H2.H1) && isscalar(H2.H1))
                this=H2;
                this.H1=this.H1*H1;
            else            
            	% strcmp is different than isa because it doesn't check all
                % superclasses as well
                if strcmp(class(H1), 'LinOpAdjoint') && isequal(H1.TLinOp,H2)
                    this.isHTH = true;
                elseif strcmp(class(H2), 'LinOpAdjoint') && isequal(H2.TLinOp,H1)
                    this.isHHt = true;
                end
            end
            if (isnumeric(H1) && isscalar(H1))
                this.name=sprintf('LinOpComposition: %s --- %s',num2str(H1),H2.name);      
            else
                this.name=sprintf('LinOpComposition: %s --- %s',H1.name,H2.name);      
            end
        end
    end
        
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
    % - applyHtH_(this,x)
    % - applyHHt_(this,x)
    % - applyAdjointInverse_(this,x)
    methods (Access = protected)
		function y = apply_(this,x) 
            % Reimplemented from :class:`LinOp`  
			if this.isHTH
				y = this.H2.applyHtH(x);
			elseif this.isHHt
				y = this.H1.applyHHt(x);
			elseif this.isH1Scalar
				y = this.H1*this.H2.apply(x);
			else
				y = this.H1.apply(this.H2.apply(x));
			end
		end
		function y = applyAdjoint_(this,x) 
            % Reimplemented from :class:`LinOp`  
			if this.isHTH || this.isHHt
				y = this.apply(x); % because self-adjoint
			elseif this.isH1Scalar
				y = this.H2.applyAdjoint(this.H1*x);
			else
				y = this.H2.applyAdjoint(this.H1.applyAdjoint(x));
			end
		end
		function y = applyHtH_(this,x)
            % Reimplemented from :class:`LinOp`  
			if this.isHTH || this.isHTH 
				y = this.apply(this.apply(x)); % because self-adjoint
			elseif this.isH1Scalar
				y = this.H2.applyAdjoint(this.H1.^2*( this.H2.apply(x)));
			else
				y = this.H2.applyAdjoint(this.H1.applyHtH( this.H2.apply(x)));
			end
        end
		function y = applyHHt_(this,x)
            % Reimplemented from :class:`LinOp`  
			if this.isHTH || this.isHTH 
				y = this.apply(this.apply(x)); % because self-adjoint
			elseif this.isH1Scalar
				y = this.H1*(this.H2.applyHHt( this.H1*x));
			else
				y = this.H1.apply(this.H2.applyHHt( this.H1.applyAdjoint(x)));
			end
        end
        function y = applyAdjointInverse_(this,x) 
            % Reimplemented from :class:`LinOp`  
            if this.isinvertible
                if this.isH1Scalar
                    y = this.H2.applyAdjointInverse(x/this.H1);
                else
                    y = this.H2.applyAdjointInverse(this.H1.applyAdjointInverse(x));
                end
            else
                x = applyAdjointInverse_@LinOp(y);
            end
        end
    end
end

