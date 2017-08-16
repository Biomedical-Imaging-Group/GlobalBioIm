classdef CostL2Composition <  CostComposition
    % CostL2Composition: Composition of a :class:`CostL2` with a
    % :class:`Map`
    % $$C(\\mathrm{x}) := \\frac12\\|\\mathrm{Hx} - \\mathrm{y}\\|^2_W $$
    %
    % :param H1: :class:`CostL2` object
    % :param H2:  :class:`Map` object
    %
    % All attributes of parent class :class:`CostL2` and :class:`CostComposition` are inherited. 
    %
    % **Example** C=CostL2Composition(H1,H2)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`CostL2`, :class:`CostComposition` :class:`LinOp`

    %%    Copyright (C) 2017 
    %     E. Soubies emmanuel.soubies@epfl.ch & 
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
    
    %% Constructor
    methods
        function this = CostL2Composition(H1,H2)
            assert(isa(H1,'CostL2'),'First argument must be a CostL2 object');
            this@CostComposition(H1,H2);
            if this.H2.norm>=0
                this.lip=this.H1.lip*this.H2.norm^2;
            end
            this.name=sprintf('CostL2Composition ( %s )',H2.name);   
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyGrad_(this,x)
    % - applyProx_(this,z,alpha)
    % - makeComposition_(this,G)
    methods (Access = protected)
        function g=applyGrad_(this,x)
        	% Reimplemented from parent class :class:`CostComposition`.
        	% $$ \\nabla C(\\mathrm{x}) = \\mathrm{J_{H}^* W (Hx - y)} $$
        	% It is L-Lipschitz continuous with \\( L \\leq \\|\\mathrm{H}\\|^2 \\|\\mathrm{W}\\|\\).
            %
            % If :attr:`doPrecomputation` is true, \\(\\mathrm{W}\\) is a scaled identity and \\(\\mathrm{H}\\)
            % is a :class:`LinOp`, then  \\(\\mathrm{H^* Wy} \\) is preconmputed and the gradient
            % computed using the :meth:`applyHtH` method.
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
        	%
            % Implemented if the operator \\(\\alpha\\mathrm{H^{\\star}WH + I}  \\) is invertible:
            % $$ \\mathrm{y} = (\\alpha\\mathrm{H^{\\star}WH + I} )^{-1} (\\alpha \\mathrm{H^TWy +x})$$
            %
            % If :attr:`doPrecomputation` is true, then
            % \\(\\mathrm{H^TWy}\\) is stored.
            if this.isH2LinOp
                if isnumeric(this.H1.W) || (isa(this.H1.W,'LinOpDiag') && this.H1.W.isScaledIdentity)
                    HtHplusId=alpha*this.H1.W*(this.H2'*this.H2) + LinOpDiag(this.H2.sizein,1);
                else
                    HtHplusId=alpha*this.H2'*this.H1.W*this.H2 + LinOpDiag(this.H2.sizein,1);
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
                end
            else
                y=applyProx_@CostComposition(this,x,alpha);
            end
% --- Old version: (in case that computing HthplusId too heavy)  
% We can maybe put the code above as default and then keep the old way 
% with the if
            % DOC
        	% - the operator :attr:`H`  is a :class:`LinOpConv` and :attr:`W`  is a :class:`LinOpIdentity`;
        	% $$\\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathcal{F}^{-1}\\left(\\frac{\\mathcal{F}(\\mathrm{x}) + \\alpha  \\mathcal{F}(\\mathrm{H}^*)\\mathcal{F}(\\mathrm{y})  }{1+\\alpha \\vert\\mathcal{F}(\\mathrm{H})\\vert^2} \\right)$$
        	% where \\(\\mathcal{F} \\) stands for the Fourier transform.
%             if isa(this.H2,'LinOpConv') && (isnumeric(this.H1.W) || isa(this.H1.W,'LinOpScaledIdentity')) 
% CODE
%                 if isnumeric(this.W), w=this.w; else, w=this.W.nu; end
%                 if this.doPrecomputation
%                     if ~isfield(this.precomputeCache,'fftHstardata')
%                         this.precomputeCache.fftHstardata=conj(this.H2.mtf).*Sfft(this.H1.y*w,this.H2.Notindex);
%                     end
%                     y=iSfft((Sfft(x,this.H2.Notindex) + w*alpha*this.precomputeCache.fftHstardata)./(1+w*alpha*(abs(this.H2.mtf).^2)), this.H2.Notindex);
%                     if ~this.H2.isComplexOut, y=real(y);end
%                 else
%                     fftHstardata=conj(this.H2.mtf).*Sfft(w*this.H1.C1y,this.H2.Notindex);
%                     y=iSfft((Sfft(x,this.H2.Notindex) + w*alpha*fftHstardata)./(1+w*alpha*(abs(this.H2.mtf).^2)), this.H2.Notindex);
%                 end
%             else
%                 y=applyProx_@CostComposition(this,x,alpha);
%             end
        end
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`Cost`
            M=CostL2Composition(this.H1,this.H2*G);
        end
    end
end
