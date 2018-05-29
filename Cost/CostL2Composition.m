classdef CostL2Composition <  CostComposition
    % CostL2Composition: Composition of a :class:`CostL2` with a
    % :class:`Map`
    % $$C(\\mathrm{x}) := \\frac12\\|\\mathrm{Hx} - \\mathrm{y}\\|^2_W = \\frac12 (\\mathrm{Hx} - \\mathrm{y})^T W (\\mathrm{Hx} - \\mathrm{y}) $$
    %
    % :param H1: :class:`CostL2` object
    % :param H2:  :class:`Map` object
    %
    % All attributes of parent class :class:`CostL2` and :class:`CostComposition` are inherited. 
    %
    % **Example** C=CostL2Composition(H1,H2)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostL2`, :class:`CostComposition`, :class:`LinOp`

    %%    Copyright (C) 2017 
    %     E. Soubies emmanuel.soubies@epfl.ch & 
    %     F. Soulez ferreol.soulez@epfl.ch &
	%     M. T. McCann michael.mccann@epfl.ch
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
    
    properties (Access=public)
        OpSumP;  % Operator used to compute the prox when H2 is the combination of a LinOpDownsample with a LinOpConv
        LLt;     % averaged convolution kernel used to compute the prox when H2 is the combination of a LinOpDownsample with a LinOpConv
    end
    
    %% Constructor
    methods
        function this = CostL2Composition(H1,H2)
            assert(isa(H1,'CostL2'),'First argument must be a CostL2 object');
            this@CostComposition(H1,H2);
            if this.H2.norm>=0
                this.lip=this.H1.lip*this.H2.norm^2;
            end
            this.name=sprintf('CostL2Composition ( %s )',H2.name);
            if isa(this.H2,'LinOpComposition') && isa(this.H2.H1,'LinOpDownsample') &&  isa(this.H2.H2,'LinOpConv') && isnumeric(this.H1.W) && this.H1.W==1
                this.OpSumP=LinOpSumPatches(this.H2.H1.sizein,this.H2.H1.sizein./this.H2.H1.df);
                this.LLt=this.OpSumP*(abs(this.H2.H2.mtf).^2);            
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this, x)
    % - applyGrad_(this,x)
    % - applyProx_(this,z,alpha)
    % - makeComposition_(this,G)
    methods (Access = protected)
		function y = apply_(this, x)
			% Reimplemented from parent class :class:`CostComposition`.
            % $$ C(\\mathrm{x}) = \\frac12\\|\\mathrm{Hx} - \\mathrm{y}\\|^2_W  $$
			%
			% If :attr:`doPrecomputation` is true, \\(\\mathrm{W}\\) is a scaled identity and \\(\\mathrm{H}\\)
			% is a :class:`LinOp`, then  \\(\\mathrm{H^* Wy} \\) and \\(\\|
			% \\mathrm{y} \\|^2_{\\mathrm{W}}\\) are
			% precomputed and  \\(C(\\mathrm{x}) \\) is evaluated using the
			% :meth:`applyHtH` method, i.e.
            % $$ C(\\mathrm{x}) = \\frac12 \\langle \\mathrm{W H^{\\star}Hx,x} \\rangle - \\langle \\mathrm{x}, \\mathrm{H^* Wy}  \\rangle  + \\frac12\\| \\mathrm{y} \\|^2_{\\mathrm{W}}$$
			if this.isH2LinOp && (isnumeric(this.H1.W) || isa(this.H1.W,'LinOpScaledIdentity')) && this.doPrecomputation
				if ~isfield(this.precomputeCache,'WHty')
					this.precomputeCache.WHty=this.H1.W*this.H2.applyAdjoint(this.H1.y);
				end
				if ~isfield(this.precomputeCache,'ytWy')
					this.precomputeCache.ytWy= this.H1.y(:)' * reshape(this.H1.W * this.H1.y, [], 1);
				end
				y=0.5 * x(:)' * reshape(this.H1.W*this.H2.applyHtH(x), [],1) - x(:)' * this.precomputeCache.WHty(:) + 0.5 * this.precomputeCache.ytWy;
				y = real(y); % due to numerical error
			else
				y=apply_@CostComposition(this,x);
			end
        end		
		function g=applyGrad_(this,x)
			% Reimplemented from parent class :class:`CostComposition`.
			% $$ \\nabla C(\\mathrm{x}) = \\mathrm{J_{H}^* W (Hx - y)} $$
            %
            % If :attr:`doPrecomputation` is true, \\(\\mathrm{W}\\) is a scaled identity and \\(\\mathrm{H}\\)
            % is a :class:`LinOp`, then  \\(\\mathrm{H^* Wy} \\) is precomputed and the gradient
            % is evaluated using the :meth:`applyHtH` method, i.e.
            % $$ \\nabla C(\\mathrm{x}) = \\mathrm{W H^{\\star}Hx} -  \\mathrm{H^* Wy}  $$
			
            if this.isH2LinOp && (isnumeric(this.H1.W) || isa(this.H1.W,'LinOpScaledIdentity')) && this.doPrecomputation
                if ~isfield(this.precomputeCache,'WHty')
                    this.precomputeCache.WHty=this.H1.W*this.H2.applyAdjoint(this.H1.y);
                end
                g=this.H1.W*this.H2.applyHtH(x) - this.precomputeCache.WHty;
            else
                g=applyGrad_@CostComposition(this,x);
            end
        end
        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`CostComposition`.
            % Implemented 
            %  - if the operator \\(\\alpha\\mathrm{H^{\\star}WH + I}  \\) is invertible:
            %    $$ \\mathrm{y} = (\\alpha\\mathrm{H^{\\star}WH + I} )^{-1} (\\alpha \\mathrm{H^TWy +x})$$
            %  - if \\(\\mathrm{H}\\) is a :class:`LinOpComposition`
            %    composing a :class:`LinOpDownsample` with a
            %    :class:`LinOpConv`. The implementation follows [1].
            %
            % **Note** If :attr:`doPrecomputation` is true, then \\(\\mathrm{H^TWy}\\) is stored.
            %
            % **References**
            %
            % [1] Zhao Ningning et al. "Fast Single Image Super-Resolution Using a New Analytical Solution for l2-l2 Problems".
            % IEEE Transactions on Image Processing, 25(8), 3683-3697 (2016).
            
            if isa(this.H2,'LinOpConv') && (isnumeric(this.H1.W) || (isa(this.H1.W,'LinOpDiag') && this.H1.W.isScaledIdentity))
            % If the composed operator is a LinOpConv
                if this.doPrecomputation
                    if ~isfield(this.precomputeCache,'fftHstardata')
                        this.precomputeCache.fftHstardata=conj(this.H2.mtf).*Sfft(this.H1.y*this.H1.W,this.H2.Notindex);
                    end
                    y=iSfft((Sfft(x,this.H2.Notindex) + this.H1.W*alpha*this.precomputeCache.fftHstardata)./(1+this.H1.W*alpha*(abs(this.H2.mtf).^2)), this.H2.Notindex);                  
                else
                    fftHstardata=conj(this.H2.mtf).*Sfft(this.H1.W*this.H1.y,this.H2.Notindex);
                    y=iSfft((Sfft(x,this.H2.Notindex) + this.H1.W*alpha*fftHstardata)./(1+this.H1.W*alpha*(abs(this.H2.mtf).^2)), this.H2.Notindex);
                end
                if this.H2.isReal, y=real(y);end
            % If the composed operator is a composition between a LinOpDownsample and a LinOpConv    
            elseif isa(this.H2,'LinOpComposition') && isa(this.H2.H1,'LinOpDownsample') &&  isa(this.H2.H2,'LinOpConv') && isnumeric(this.H1.W) && this.H1.W==1
                 % this.H1 -> CostL2
                 % this.H2.H1 -> LinOpDownsample
                 % this.H2.H2 -> LinOpConv
                 if this.doPrecomputation
                    if ~isfield(this.precomputeCache,'fftHstarSdata')
                        this.precomputeCache.fftHstarSdata=conj(this.H2.H2.mtf).*Sfft(this.H2.H1.applyAdjoint(this.H1.y),this.H2.H2.Notindex);
                    end
                    fftr=this.precomputeCache.fftHstarSdata+Sfft(x/alpha,this.H2.H2.Notindex);
                    y=real(iSfft(alpha*(fftr - conj(this.H2.H2.mtf).*this.OpSumP.applyAdjoint(this.OpSumP.apply(this.H2.H2.mtf.*fftr)./(prod(this.H2.H1.df)/alpha+this.LLt))),this.H2.H2.Notindex));
                 else
                    fftr=conj(this.H2.H2.mtf).*Sfft(this.H2.H1.applyAdjoint(this.H1.y),this.H2.H2.Notindex)+Sfft(x/alpha,this.H2.H2.Notindex);
                    y=real(iSfft(alpha*(fftr - conj(this.H2.H2.mtf).*this.OpSumP.applyAdjoint(this.OpSumP.apply(this.H2.H2.mtf.*fftr)./(prod(this.H2.H1.df)/alpha+this.LLt))),this.H2.H2.Notindex));
                 end                                           
            % Default implementation
            elseif this.isH2LinOp && ~this.isH2SemiOrtho
                if ~this.doPrecomputation || (this.doPrecomputation && (~isfield(this.precomputeCache,'HtHplusId')  || (alpha~=this.precomputeCache.alpha)))
                    if isnumeric(this.H1.W) || (isa(this.H1.W,'LinOpDiag') && this.H1.W.isScaledIdentity)
                        HtHplusId=alpha*this.H1.W*(this.H2'*this.H2) + LinOpDiag(this.H2.sizein,1);
                    else
                        HtHplusId=alpha*this.H2'*this.H1.W*this.H2 + LinOpDiag(this.H2.sizein,1);
                    end
                    if this.doPrecomputation
                        this.precomputeCache.HtHplusId=HtHplusId;
                        this.precomputeCache.alpha=alpha;
                    end
                else
                    HtHplusId=this.precomputeCache.HtHplusId;
                end
                if HtHplusId.isInvertible
                    if this.doPrecomputation
                        if ~isfield(this.precomputeCache,'HtWy')
                            this.precomputeCache.HtWy=this.H2.applyAdjoint(this.H1.W*this.H1.y);
                        end
                        y=HtHplusId.applyInverse(alpha*this.precomputeCache.HtWy+x);
                    else
                        y=HtHplusId.applyInverse(alpha*this.H2.applyAdjoint(this.H1.W*this.H1.y)+x);
                    end
                else
                    y=applyProx_@CostComposition(this,x,alpha);
                end
            else
                y=applyProx_@CostComposition(this,x,alpha);
            end
            % Emmanuel: Actually (for deconv) when the doPrecomputation is
            % activated, the first if is not faster (if alpha do not change)
        end
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`Cost`. Instantiates a new
            % :class:`CostL2Composition` with the updated composed
            % :class:`Map`.
            M=CostL2Composition(this.H1,this.H2*G);
        end
    end
end
