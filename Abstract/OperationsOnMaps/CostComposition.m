classdef CostComposition < MapComposition & Cost
    % CostComposition: Compose a :class:`Cost` with a :class:`Map`
    % $$ C(\\mathrm{x}) := F( \\mathrm{Hx}) $$
    % where \\(F\\) is a :class:`Cost` and \\(\\mathrm{H}\\) a
    % :class:`Map`
    %
    % :param H1:  :class:`Cost` 
    % :param H2:  :class:`LinOp`
    %
    % **Example** C = CostComposition(H1,H2)
    %
    % See also :class:`Map`, :class:`MapComposition`, :class:`Cost`

    %%    Copyright (C) 2017 
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
 
    properties  (SetAccess = protected,GetAccess = protected)
        isH2LinOp=false;
        isH2SemiOrtho=false;
        nu=0;
    end
    properties
        isInputReal=0;% to correctly switch in the  grad
    end
    
    %% Constructor
    methods
        function this = CostComposition(H1,H2)
            this@MapComposition(H1,H2);   
            assert(isa(H1,'Cost'),'First argument should be a Cost');
            if isa(H2,'LinOp') && H1.isConvex
                this.isConvex=true;
            end 
            if isa(H2,'LinOp')
                this.isH2LinOp=true;
                T=this.H2*this.H2';
                if isa(T,'LinOpDiag') && T.isScaledIdentity
                    if T.diag>0
                        this.isH2SemiOrtho=true;
                        this.nu=T.diag;
                    end
                end
                this.isSeparable=(H1.isSeparable && (isa(H2,'LinOpDiag') || isa(H2,'LinOpSelector')));
            end
            this.name=sprintf('CostComposition( %s ; %s )',H1.name,H2.name);  
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyGrad_(this,x)
    % - applyProx_(this,z,alpha)
    % - makeComposition_(this,G)
    methods (Access = protected)
        function g=applyGrad_(this,x)
            % Reimplemented from :class:`Cost`
            g = this.H2.applyJacobianT(this.H1.applyGrad(this.H2.apply(x)),x);
            if this.isInputReal
                g = real(g);
            end
        end
        function x=applyProx_(this,z,alpha)
            % Reimplemented from :class:`Cost`
            %
            % If this.H2 is a :class:`LinOp` and \\(\\mathrm{H}
            % \\mathrm{H}^{\\star}\\) is a :class:`LinOpScaledIdentity`
            % 
            if this.isConvex && this.isH2LinOp && this.isH2SemiOrtho 
                H2z=this.H2*z;
                x = z + 1/this.nu*this.H2.applyAdjoint(this.H1.applyProx(H2z,alpha*this.nu)-H2z);
            else
                x = applyProx_@Cost(this,z,alpha);
            end
        end
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`Cost`
            if isa(G,'LinOp')
                T=G*G';
                if isa(T,'LinOpDiag') && T.isScaledIdentity
                    M=CostComposition(this,G);
                else
                    M=CostComposition(this.H1,this.H2*G);
                end
            else
                M=CostComposition(this.H1,this.H2*G);
            end
        end
    end
    
    %% Utility methods
    % - set_y(this,y)
    methods
        function set_y(this,y)
            % Set the attribute \\(\\mathrm{y}\\)
            %
            %  - has to be conformable with the :attr:`sizeout` of the
            %    :class:`Map`\\(\\mathrm{H}\\),
            %  - can be anything if \\(\\mathrm{H}\\) is not yet set (empty),
            %  - can be a scalar.
            assert(isnumeric(y),'y must be a numeric');
            assert(isscalar(y) || checkSize(y,this.H2.sizeout),'size y must be a scalar or equal to this.H2.sizeout');
            this.y=y;
        end
    end
    
    methods (Access = protected)
         %% Copy
         function this = copyElement(obj)
             this = copyElement@MapComposition(obj);
            this.memoCache.applyGrad=struct('in', [], 'out', []);
            this.memoCache.applyProx=struct('in', [], 'out', []);
            this.memoCache.applyProxFench=struct('in', [], 'out', []);
      end
    end  
end
