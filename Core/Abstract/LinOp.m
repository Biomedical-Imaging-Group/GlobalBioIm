classdef (Abstract) LinOp < Map
    % Abstract class for linear operators
	% $$ \\mathrm{H}: \\mathrm{X}\\rightarrow \\mathrm{Y}.$$
    %
    % :param name: name of the linear operator \\(\\mathbf{H}\\)
    % :param sizein:  dimension of the left hand side vector space \\(\\mathrm{X}\\) 
    % :param sizeout:  dimension of the right hand side vector space \\(\\mathrm{Y}\\) 
    % :param isInvertible: true if the operator is invertible
    % :param isComplex: true if the operator is complex
    % :param norm: norm of the operator \\(\\|\\mathrm{H}\\|\\) (if known, otherwise -1)
	%
 
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
    
	properties
		
		%memoizeOpts = struct('apply', false,... % todo: can't change the default value in this way
		%	'adjoint', false); % todo etc
	end
	properties (SetAccess = private)
			%memoCache = struct('apply', struct('in', [], 'out', [])); % etc
	end
	
	methods (Sealed)
		function x = applyAdjoint(this, y)
			% x = H.applyAdjoint(y) returns the application of the adjoint of
			% the LinOp H to the vector y. 
			%
			
			% check input size
			if ~checkSize(y, this.sizeout)
				error('Input to applyAdjoint was size [%s], didn''t match stated sizeout: [%s].',...
					num2str(size(y)), num2str(this.sizeout));
			end
			
			% memoize
			x = this.memoize('applyAdjoint', @this.applyAdjoint_, y);
			
			% check output size
			if ~checkSize(x, this.sizein)
				warning('Output of applyAdjoint was size [%s], didn''t match stated sizein: [%s].',...
					num2str(size(x)), num2str(this.sizein));
			end
			
		end
		
	  
	  function this = transpose(this)
		% Overload operator (.') for :class:`LinOp` objects (i.e. :meth:`adjoint`).
		%
		% **Note**: operators (.') and (') (see :meth:`ctranspose`) are identical for :class:`LinOp`.
		
		if this.isComplex
		  warning('Warning: Do you mean adjoint? For LinOp object transpose() is an alias of adjoint method');
		end
		this = Adjoint(this);
	  end
	  
	  function this = ctranspose(this)
		% Overload operator (') for :class:`LinOp` objects (i.e. :meth:`adjoint`).
		%
		% **Note**: operators (') and (.') (see :meth:`transpose`) are identical for :class:`LinOp`.
		
		this = Adjoint(this);
	  end
	  
	  function Hnew = mtimes(this,H2)
		% Overload operator (*) for :class:`LinOp` objects.
		% $$ \\mathrm{H}_{new} =  \\mathrm{H_2} \\mathrm{H}$$
		%
		% :param H2: :class:`LinOp` object or a scalar in \\(\\mathbb{R}\\)
		% :returns Hnew: :class:`LinOp`
		
		if isa(H2,'LinOp')
		  Hnew = MulLinOp(this, H2);
		else
		  Hnew = this.apply(H2);
		end
	  end
	  
	  function Hnew = plus(this,H2)
		% Overload operator (+) for :class:`LinOp` objects.
		% $$ \\mathrm{H}_{new} =  \\mathrm{H_2} + \\mathrm{H}$$
		%
		% :param H2: :class:`LinOp` object
		% :returns Hnew: :class:`LinOp`
		
		assert(isa(H2,'LinOp'),'addition of LinOp is only define with other LinOp');
		Hnew = SumLinOp({this,H2});
	  end
	  
	  
	  
	  %HtH
	  
	  %HHt
	end
	
	methods (Access = protected, Sealed)
	  function x = applyJacobianT_(this, y, v)
		x = this.applyAdjoint_(y);
	  end
	end
	  
    methods (Access = protected) % all the other underscore methods
	  
        function x = adjoint_(this, y) 
            % **(Abstract method)** Apply the adjoint of the linear operator
        	%
        	% :param x: \\(\\in Y\\)
        	% :returns y: \\(= \\mathrm{H^*x}\\) where \\(\\mathrm{H}^*\\) is defined such that $$\\langle\\mathrm{Hx},\\mathrm{y}\\rangle_{\\mathrm{Y}} = \\langle \\mathrm{x}, \\mathrm{H}^*\\mathrm{y} \\rangle_{\\mathrm{X}}$$
        	
			% todo: general case possible
		end
		
		function x = applyAdjoint_(this, y)
			% **(Abstract method)** Apply the adjoint of the linear operator
			%
			% :param x: \\(\\in Y\\)
			% :returns y: \\(= \\mathrm{H^*x}\\) where \\(\\mathrm{H}^*\\) is defined such that $$\\langle\\mathrm{Hx},\\mathrm{y}\\rangle_{\\mathrm{Y}} = \\langle \\mathrm{x}, \\mathrm{H}^*\\mathrm{y} \\rangle_{\\mathrm{X}}$$
			
			error('adjoint not implemented');
		end
		
        
        function y = HtH_(this,x) 
            % Apply \\(\\mathrm{H}^*\\mathrm{H}\\)
        	%
        	% :param x: \\(\\in X\\)
        	% :returns y: \\(= \\mathrm{H^*Hx}\\)
        	%
        	% **Note**: There is a default implementation in the abstract class :class:`LinOp` which calls successively the :meth:`apply` and :meth:`adjoint` methods. However, it can be reimplemented in derived classes if there exists a faster way to perform computation.
        	
            y = this.applyAdjoint(this.apply(x));
        end
        
        function y = HHt_(this,x) 
            % Apply \\(\\mathrm{H}\\mathrm{H}^*\\)
        	%
        	% :param x: \\(\\in Y\\)
        	% :returns y: \\(= \\mathrm{HH^*x}\\)
        	%
        	% **Note**: There is a default implementation in the abstract class :class:`LinOp` which calls successively the :meth:`adjoint` and :meth:`apply` methods. However, it can be reimplemented in derived classes if there exists a faster way to perform computation.
        	
            y = this.apply(this.applyAdjoint(x));
		end
		
		function G = makeComposition_(this, H)
			if isa(H, 'LinOp')
				G = MulLinOp({this, H});
			else
				G = makeComposition_@Map(this, H);
			end
				
		end
        
		% todo: ask ferreol: is this needed? when do we have fast HT W H?
%         function y = HtWH(this,x,W) 
%         	% Apply \\(\\mathrm{H}^*\\mathrm{WH}\\)
%         	%
%         	% :param x: \\(\\in X\\)
%         	% :param W: a :class:`LinOp` object
%         	% :returns y: \\(= \\mathrm{H^*WHx}\\)
%         	%
% 
%             if (isscalar(W) && isreal(W))
%                 y = W.*this.HtH(x);
%             else
%                 assert(isa(W,'LinOp'),'W must be a LinOp');
%                 y = this.adjoint(W.apply(this.apply(x)));
%             end
%         end
        
      
        
        function y= adjointInverse(this,x) 
            % **(Abstract method)** Apply \\(\\mathrm{H}^{-*}\\) (if applicable)
        	%
        	% :param x: \\(\\in X\\)
        	% :returns y: \\(= \\mathrm{H^{-*}x}\\)
        	%
        	
            if this.isInvertible
                error('adjointInverse not implemented');
            else
                error('Operator not invertible');
            end
		end
        
	end
		
		methods 
		
   
		
	end
		

end

