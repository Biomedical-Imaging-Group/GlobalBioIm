classdef LinOpAdjoint < LinOp
    % Adjoint : Builds the adjoint LinOp
    %
    % :param TLinOp: :class:`LinOp` object
    %
    % **Example** Tadj=LinOpAdjoint(TLinOp)
    %
    % See also :class:`Map`, :class:`LinOp`
    
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
        TLinOp     % linop
    end
    
    %% Constructor
    methods 
        function this = LinOpAdjoint(TLinOp)
            this.name =sprintf('%s Adjoint', TLinOp.name);           
            assert(isa(TLinOp,'LinOp'),'Input should be a  LinOp');
          
            this.TLinOp = TLinOp;
            this.isDifferentiable= this.TLinOp.isDifferentiable;
            this.isInvertible=this.TLinOp.isInvertible;
            this.sizein =  this.TLinOp.sizeout;
            this.sizeout =  this.TLinOp.sizein;			
			this.norm = this.TLinOp.norm;  
             
        end       
    end
	
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
    % - applyHtH_(this,x)
    % - applyHHt_(this,x)
    % - applyInverse_(this,x)
    % - applyAdjointInverse_(this,x)
    % - makeAdjoint_(this)
    % - makeHHt_(this)
    % - makeHtH_(this)
	methods (Access = protected)
        function y = apply_(this,x) 
            % Reimplemented from :class:`LinOp`
            y =this.TLinOp.applyAdjoint(x);
        end
        function y = applyAdjoint_(this,x) 
            % Reimplemented from :class:`LinOp`
            y =this.TLinOp.apply(x);
        end
        function y = applyHtH_(this,x)
            % Reimplemented from :class:`LinOp`
            y =this.TLinOp.HHt(x);
        end
        function y = applyHHt_(this,x)
            % Reimplemented from :class:`LinOp`
            y =this.TLinOp.HtH(x);
        end
        function y = applyInverse_(this,x)
            % Reimplemented from :class:`LinOp`
            y =this.TLinOp.applyAdjointInverse(x);
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from :class:`LinOp`
            y =this.TLinOp.applyInverse(x); 
        end
        function M = makeAdjoint_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=this.TLinOp;
        end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=this.TLinOp.makeHtH();
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=this.TLinOp.makeHHt();
        end
    end
    
    methods (Access = protected)
        %% Copy
        function this = copyElement(obj)
            this = copyElement@LinOp(obj);
            this.TLinOp = copyElement@LinOp(obj.TLinOp);
        end
    end
end

