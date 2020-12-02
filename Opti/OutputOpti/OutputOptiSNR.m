classdef OutputOptiSNR  < OutputOpti
    % OutputOptiSNR class for algorithms displayings and savings
    %
    % At each :attr:`ItUpOut` iterations of an optimization algorithm (see :class:`Opti` generic class),
    % the update method of an :class:`OutputOptiSNR` object will be executed in order to acheive user
    % defined computations, e.g.,
    %
    %  - compute cost / SNR
    %  - store current iterate / cost value
    %  - plot/display stuffs
    %
    % The present generic class implements a basic update method that:
    %
    %  - display the iteration number
    %  - computes & display thhandlee cost (if activated)
    %  - computes & display the SNR if ground truth is provided
    %
    % :param name:  name of the :class:`OutputOptiSNR`
    % :param computecost:  boolean, if true the cost function will be computed
    % :param xtrue: ground truth to compute the error with the solution (if provided)
    % :param evolcost: array to save the evolution of the cost function
    % :param evolsnr: array to save the evolution of the SNR
    % :param saveXopt: boolean (defaul true) to save the evolution of the optimized variable xopt.
    % :param evolxopt:  cell saving the optimization variable xopt
    % :param iterVerb:  message will be displayed every iterVerb iterations (must be a multiple of the :attr:`ItUpOut` parameter of classes :class:`Opti`)
    % :param costIndex: select a specific cost function among a sum in the case where the optimized cost function is a sum of cost functions
    %
    % **Example** OutOpti=OutputOptiSNR(computecost,xtrue,iterVerb)
    %
    % **Important** The update method should have an unique imput that is the :class:`Opti` object in order to
    % be generic for all Optimization routines. Hence the update method has access (in reading mode)
    % to all the properties of :class:`Opti` objects.
    %
    % See also :class:`Opti`
    
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
    
    properties (SetAccess = protected,GetAccess = public)
        evolsnr;           % array saving the evolution of the error with the groud truth
        xtrue;             % Ground Truth
        snrOp=[];          % Map to which the coefficients must be applied to compute reconstructed signal
        normXtrue;         % Norm of the true signal (to compute snr)
    end
    
    methods
        %% Constructor
        function this=OutputOptiSNR(computecost,xtrue,iterVerb,costIndex,snrOp)
            
            this.name = 'OutputOptiSNR';% name
            if nargin>=1
                if isscalar(computecost)
                    computecost = (computecost ~= 0);
                end
                
                assert(islogical(computecost),'Parameter computecost must be logical');
                this.computecost=computecost;
                
            end
            if nargin>=2
                this.xtrue=xtrue;
            end
            if nargin>=3
                assert(isscalar(iterVerb) && iterVerb>=0,'Parameter iterVerb must be a positive integer');
                this.iterVerb=iterVerb;
            end
            if nargin>=4, this.costIndex=costIndex;end
            if ~isempty(this.xtrue)
                this.xtrue=xtrue;
                this.normXtrue=norm(this.xtrue(:));
            end
            if nargin==5, this.snrOp = snrOp;end
            if ~isempty(this.snrOp)
                assert(~isempty(this.xtrue), 'Ground truth must be provided in order to compute SNR');
                assert(isa(this.snrOp, 'Map') && isequal(this.snrOp.sizeout, size(this.xtrue)), ...
                    'The SNR operator must have the same output size as the ground truth');
            elseif ~isempty(this.xtrue)
                this.snrOp = LinOpDiag(size(this.xtrue)); % Identity operator by default
            end
        end
        %% Initialization
        function init(this)
            % Initialize the arrays and counters.
            this.count=1;
            this.evolcost=zeros_(1);
            this.evolsnr=zeros_(1);
            this.iternum = [];
            this.evolxopt = {};
        end
        %% Update method
        function update(this,opti)
            % Computes SNR, cost and display evolution.
            str=sprintf('Iter: %5i',opti.niter);
            if this.computecost
                cc=this.computeCost(opti);
                str=sprintf('%s | Cost: %4.4e',str,cc);
                this.evolcost(this.count)=cc;
            end
            
            snr=this.computeSNR(opti);
            str=sprintf('%s | SNR: %4.4e dB',str,snr);
            this.evolsnr(this.count)=snr;
            
            if this.saveXopt
                this.evolxopt{this.count}=opti.xopt;
            end
            this.iternum(this.count)=opti.niter;
            this.count=this.count+1;
            if opti.verbose && (opti.niter~=0 && (mod(opti.niter,this.iterVerb)==0) || (opti.niter==1 && this.iterVerb~=0))
                disp(str);
            end
        end
        function cc=computeCost(this,opti)
            % Evaluate the cost function at the current iterate xopt of
            % the given :class:`Opti` opti object
            if (any(this.costIndex>0) && isa(opti.cost,'CostSummation'))
                cc = 0;
                for n=1:numel(this.costIndex)
                    cc = cc+opti.cost.mapsCell{this.costIndex(n)}*opti.xopt;
                end
            else
                cc = opti.cost*opti.xopt;
            end
        end
        function snr=computeSNR(this,opti)
            % Evaluate the snr for the current iterate xopt of
            % the given :class:`Opti` opti object
            reconstruction = this.snrOp.apply(opti.xopt);
            snr=20*log10(this.normXtrue/norm(this.xtrue(:)-reconstruction(:)));
        end
    end
end
