classdef LinOpComposition <  LinOp & MapComposition
    % LinOpComposition : Composition of LinOps
    % $$ \\mathrm{H}(\\mathrm{x}) = \\mathrm{H}_1 \\mathrm{H}_2\\mathrm{x} $$
    %
    % :param H1:  left hand side :class:`LinOp` (or a scalar)
    % :param H2:  right hand side :class:`LinOp`
    %
    % **Example** H=LinOpComposition(H1,H2)
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
    
    %% Properties
    % - Private
    properties (SetAccess = protected,GetAccess = protected)
		isHTH = 0;
		isHHt = 0;
    end
    
    %% Constructor 
    methods
        function this = LinOpComposition(H1,H2)
            % Call superclass constructor
            this@MapComposition(H1,H2); 
            % Initialize
            this.initObject('LinOpComposition');
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclasses :class:`MapComposition` and :class:`LinOp`
            
            % Call superclass method 
            updateProp@MapComposition(this,prop);
            updateProp@LinOp(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'H1')  || strcmp(prop,'H2') ||  strcmp(prop,'all')
                % strcmp is different than isa because it doesn't check all
                % superclasses as well
                assert(isa(this.H2,'LinOp') && isa(this.H1,'LinOp'),'H1 and H2 have to be a LinOps');
                if strcmp(class(this.H1), 'LinOpAdjoint') && isequal(this.H1.TLinOp,this.H2)
                    this.isHTH = true;
                elseif strcmp(class(this.H2), 'LinOpAdjoint') && isequal(this.H2.TLinOp,this.H1)
                    this.isHHt = true;
                end
                this.name=sprintf('LinOpComposition( %s ; %s )',this.H1.name,this.H2.name);
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
            else
                y = this.H1.apply(this.H2.apply(x));
            end
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from :class:`LinOp`
            if this.isHTH || this.isHHt
                y = this.apply(x); % because self-adjoint
            else
                y = this.H2.applyAdjoint(this.H1.applyAdjoint(x));
            end
        end
        function y = applyHtH_(this,x)
            % Reimplemented from :class:`LinOp`
            if this.isHTH || this.isHTH
                y = this.apply(this.apply(x)); % because self-adjoint
            else
                y = this.H2.applyAdjoint(this.H1.applyHtH( this.H2.apply(x)));
            end
        end
        function y = applyHHt_(this,x)
            % Reimplemented from :class:`LinOp`
            if this.isHTH || this.isHTH
                y = this.apply(this.apply(x)); % because self-adjoint
            else
                y = this.H1.apply(this.H2.applyHHt( this.H1.applyAdjoint(x)));
            end
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from :class:`LinOp`
            if this.isinvertible
                y = this.H2.applyAdjointInverse(this.H1.applyAdjointInverse(x));
            else
                y = applyAdjointInverse_@LinOp(x);
            end
        end
        function M = makeAdjoint_(this) 
            % Reimplemented from parent class :class:`LinOp`.
                M=this.H2'* this.H1';      
        end
        function M = makeHtH_(this)
            % Reimplemented from :class:`LinOp`
            M=this.H2'*this.H1.makeHtH()*this.H2;
        end
        function M = makeHHt_(this)
            % Reimplemented from :class:`LinOp`
            M=this.H1*this.H2.makeHHt()*this.H1';
        end
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`MapComposition`
             if ~isa(this.H1, 'LinOpComposition')  &&  isa(G, 'LinOpComposition')                 
                M =  this.H1*(this.H2*G); % Since H1*H2 is not simplified,
             elseif isa(G,'LinOp')
                 H2G = this.H2*G; % Since H1*H2 is not simplified, first try to compose H2 with G
                 if ~isa(H2G, 'LinOpComposition')
                     M=this.H1*H2G;
                 else
                     M = LinOpComposition(this, G);
                 end
                 H2G.clearListenerList({this.H2,G});
             else
                 M=makeComposition_@MapComposition(this,G);
            end
        end
    end
end

