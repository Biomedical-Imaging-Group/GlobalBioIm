classdef OptiRichLucy < Opti
    % Richardson-Lucy algorithm [1,2] which minimizes the KullbackLeibler
    % divergence :class:`CostKullLeib` (with TV regularization [3]).
    % $$ C(\\mathrm{x})= F(\\mathrm{x}) + \\lambda \\Vert \\mathrm{x} \\Vert_{\\mathrm{TV}} $$
    %
    % :param F: :class:`CostKullLeib` object or a :class:`CostComposition`
    %  with a :class:`CostKullLeib` and a :class:`LinOp`
    % :param TV: boolean true if TV regularization used  (default false)
    % :param lambda: regularization parameter (when TV used)
    % :param epsl: smoothing parameter to make TV differentiable at 0 (default \\(10^{-6}\\))
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
    % **Note** An interesting property of this algorithm is that it ensures 
    % the positivity of the solution from any positive initialization.
    % However, when TV is used, the positivity of the iterates is not ensured 
    % anymore if \\(\\lambda \\) is too large. Hence, \\(\\lambda \\) needs to be carefully chosen.
    %
    % **References**
    %
	% [1] Lucy, Leon B. "An iterative technique for the rectification of observed distributions" The astronomical journal (1974)
    %
	% [2] Richardson, William Hadley. "Bayesian-based iterative method of image restoration." JOSA (1972): 55-59.
    %
    % [3] N. Dey et al. "Richardson-Lucy Algorithm With Total Variation Regularization for 3D Confocal Microscope 
    % Deconvolution." Microscopy research and technique (2006).
    %
    % **Example** RL=OptiRichLucy(F,TV,lamb)
    %
    % See also :class:`Opti`, :class:`OutputOpti`, :class:`Cost`,
    % :class:`CostKullLeib`

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
    
    % Full public properties     
    properties (SetAccess = public,GetAccess = public)
		epsl=1e-6; % smoothing parameter for TV
    end
    % Protected set public read    
    properties (SetAccess = protected,GetAccess = public)
		TV=false;   % boolean true if RL-TV version activated (default false)
		lamb=1e-2;  % regularization parameter
	end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
		G;   % gradient operator (used if TV activated)
		F;   % Kullback-Leibler divergence
        isH; % boolean true if F is a CostComposition
        data;
        He1;
        bet;
    end
   
    methods
    	%% Constructor
    	function this=OptiRichLucy(F,TV,lamb)
    		this.name='Opti Richardson-Lucy';
    		assert((isa(F,'CostComposition') && isa(F.H1,'CostKullLeib') && isa(F.H2,'LinOp') ) || ...
                isa(F,'CostKullLeib'), 'The minimized functional should be the FuncKullLeib or a CostComposition between a CostKullLeib and a LinOp');
            if ~isa(F,'CostComposition')
                this.F=F*LinOpDiag(F.sizein);
            else
                this.F=F;
            end
            if nargin>=2 && ~isempty(TV), this.TV=TV; end
            this.cost=this.F;
            if nargin>=3 && ~isempty(lamb), this.lamb=lamb; end
            if this.TV
                this.G=LinOpGrad(this.F.sizein);
                if length(this.F.sizein)==2  % 2D
                    this.cost=this.cost + this.lamb*CostMixNorm21([this.F.sizein,2],3)*this.G;
                elseif length(this.F.sizein)==3
                    this.cost=this.cost + this.lamb*CostMixNorm21([this.F.sizein,3],4)*this.G;
                end
            end
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if ~any(this.xopt(:)), this.xopt=ones(size(this.xopt)); end
            this.data=this.F.H1.y;
            this.bet=this.F.H1.bet;
            this.He1=max(this.F.H2.applyAdjoint(ones_(this.F.H2.sizeout)),this.bet);
            if this.bet==0, error('Smoothing parameter beta has to be different from 0 (see constructor of CostKullLeib)'); end;
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`. For details see [1-3].
            
            if ~this.TV
                this.xopt=this.xopt./this.He1.*this.F.H2.applyAdjoint(this.data./(this.F.H2.apply(this.xopt)+this.bet));
            else
                tmp=this.G.apply(this.xopt);
                if length(size(tmp))==2     % 1D
                    nor=sqrt(tmp.^2+this.epsl);
                elseif length(size(tmp))==3 % 2D
                    nor=repmat(sqrt(sum(tmp.^2,3)+this.epsl),[1,1,size(tmp,3)]);
                elseif length(size(tmp))==4 % 3D
                    nor=repmat(sqrt(sum(tmp.^2,4)+this.epsl),[1,1,1,size(tmp,4)]);
                end
                gradReg=this.G.applyAdjoint(tmp./nor);
                this.xopt=this.xopt./(this.He1 + this.lamb*gradReg).*this.F.H2.applyAdjoint(this.data./(this.F.H2.apply(this.xopt)+this.bet));
                if sum(this.xopt(:)<0)~=0
                    warnStruct = warning('off','backtrace');
                    warning('Violation of the positivity of the solution (the regularization parameter should be decreased).');
                    warning(warnStruct);
                end
            end
            flag=this.OPTI_NEXT_IT;
        end
 
	end
end
