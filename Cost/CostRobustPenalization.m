classdef CostRobustPenalization < Cost
    % CostRobustPenalization: Robust penalization cost function
    % $$C(\\mathrm{x}) := sum rho((A.x-y)/s)$$ 
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param M: MAP to apply on the input (default Identity)
    % :param y: data vector (default 0)
    % :param options: structure containing the different options of the robust penalization
    %
    % :param options.method: objective function method to compute the cost function
    %   - 'Andrews'
    %   - 'Beaton-Tukey'
    %   - 'Cauchy' (default)
    %   - 'Fair'
    %   - 'Huber'
    %   - 'Logistic'
    %   - 'Talwar-Hinich'
    %   - 'Welsch-Dennis'
    %  
    % :param options.mu: parameters to tune the method. Default:
    %   - 'Andrews'       -> 1.339
    %   - 'Beaton-Tukey'  -> 4.685
    %   - 'Cauchy'        -> 2.385
    %   - 'Fair'          -> 1.400
    %   - 'Huber'         -> 1.345
    %   - 'Logistic'      -> 1.205
    %   - 'Talwar-Hinich' -> 2.795
    %   - 'Welsch-Dennis' -> 2.985
    %
    % :param options.flag_s: method to scale the residue before applying the objective function (default: none)
    %   - 'none'	-> no scaling applied
    %   - 'MAD'   -> median absolute deviation
    %   - 'STD'   -> standard deviation
    %   - numeric -> scaling (can be a matrix to have a different scaling for each variables)
    %
    % :param options.noise_model: model of the noise to scale the residues according to the input variable x
    %   - 'Poisson'   -> scaling_factor = \sqrt(x+var_0)
    %   - 'none'      -> no scaling according to a noise model
    %
    % :param options.var_0: value of the variance of the data at null flux (default 0)
    % :param options.eta: ratio to scale the model in the corresponding unit for the Poisson noise
    % :param options.flag_memoizeRes: memoize the computation of the residues? (default: true)
    %
    % **Example** C=CostRobustPenalization(M, y)
    %
    % **Example** C=CostRobustPenalization(M, y, options)
    %
    % **Example** C=CostRobustPenalization(M, y, [])
    %
    % See also :class:`Map`, :class:`Cost`
    
    %%    Copyright (C) 2018
    %     Created: 07/25/2018 (mm/dd/yyyy)
    %     Modified: 10/12/2018 (mm/dd/yyyy) Inspired from the deprecated
    %     CostReweightedL2
    %     Anthony Berdeu (Laboratoire Hubert Curien)
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
 
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        M;              % MAP to apply on the input
        method;         % Method tom compute the objective function
        mu;             % Parameters of the method
        flag_s;         % Method to compute the residue scaling
        noise_model;    % noise_model
            % 0 -> no noise model
            % 1 -> Poisson noise
        var_0;          % variance at null flux for the Poisson noise
        eta;            % ratio to the convert the model in the correct
            % for the Poisson noise
    end
    
    %% Constructor
    methods     	
        function this = CostRobustPenalization(M, y, options)
            
            %% Declaration
            this.name='CostRobustPenalization';
            assert(isa(M,'Map'),'M must be Map');
            this.M = M ;
            this.sizein = this.M.sizein ;
            this.isDifferentiable = true;
            
            %% Checking inputs
            % y
            if isempty(y)
                this.y = 0 ;
            else
                this.y = y ;
            end
            
            % options
            if nargin<3 || isempty(options) % Default values
                options.method = 'Cauchy' ;
                options.mu = [] ;
                options.flag_s = 'none' ;
                options.noise_model = 'none' ;
                options.flag_memoizeRes = true ;
            end
            
            %% Declaring cost's fields
            % method
            if isfield(options, 'method')
                this.method = options.method ;
            else
                this.method = 'Cauchy' ;
            end
            switch this.method
                case 'Andrews'
                    this.mu = 1.339 ;
                    this.isConvex = false ;
                case 'Beaton-Tukey'
                    this.mu = 4.685 ;
                    this.isConvex = false ;
                case 'Cauchy'
                    this.mu = 2.385 ;
                    this.isConvex = false ;
                case 'Fair'
                    this.mu = 1.400 ;
                    this.isConvex = true ;
                case 'Huber'
                    this.mu = 1.345 ;
                    this.isConvex = true ;
                case 'Logistic'
                    this.mu = 1.205 ;
                    this.isConvex = true ;
                case 'Talwar-Hinich'
                    this.mu = 2.795 ;
                    this.isConvex = false ;
                case 'Welsch-Dennis'
                    this.mu = 2.985 ;
                    this.isConvex = false ;
                otherwise
                    error(['Unknown method... ', ...
                        'Please choose a method among: ', ...
                        'Andrews', ...
                        ', Beaton-Tukey', ...
                        ', Cauchy', ...
                        ', Fair', ...
                        ', Huber', ...
                        ', Logistic', ...
                        ', Talwar-Hinich', ...
                        ', Welsch-Dennis', ...
                        ]) ;
            end
            
            % mu
            if isfield(options, 'mu')
                this.mu = options.mu ;
            end
            
            % flag_s
            if isfield(options,'flag_s')
                this.flag_s = options.flag_s ;
            else
                this.flag_s = 'none' ;
            end
            if ~isnumeric(this.flag_s)
                switch this.flag_s
                    case 'none'
                        this.flag_s = 1 ;
                    case 'MAD'
                        % median absolute deviation
                    case 'STD'
                        % standard deviation
                    otherwise
                        error(['Unknown method for the residues', ...
                            ' scaling... ', ...
                            'Please choose a method among: ', ...
                            'STD', ...
                            'MAD', ...
                            ', none', ...
                            ', a numerical value', ...
                            ]) ;
                end
            end
            
            % noise_model
            if isfield(options, 'noise_model')
                this.noise_model = options.noise_model ;
            else
                this.noise_model = 'none' ;
            end
            switch this.noise_model
                case 'none'
                    this.noise_model = 0 ;
                case 'Poisson'
                    this.noise_model = 1 ;
                    % var_0
                    if isfield(options, 'var_0')
                        this.var_0 = options.var_0 ;
                    else
                        warning(['The variance at null model is set ', ...
                            'to 0... /!\ This can lead to divergences', ...
                            ' of the residues...']) ;
                        this.var_0 = 0 ;
                    end
                    if isfield(options, 'eta')
                        this.eta = options.eta ;
                    else
                        warning(['The ratio to convert the model is', ...
                            ' set to 1... ']) ;
                        this.eta = 1 ;
                    end
                otherwise
                    error(['Unknown model for the noise... ', ...
                        'Please choose a model among: ', ...
                        'none, ', ...
                        'Poisson', ...
                        ]) ;
            end
            
            
            %% Memoize options
            if ~isfield(options, 'flag_memoizeRes')
                options.flag_memoizeRes = true ;
            end
            this.memoizeOpts.computeRes = options.flag_memoizeRes ;
            if this.memoizeOpts.computeRes
                this.memoCache.computeRes = struct('in', [], 'out', []);
            end
        end   
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    methods (Access = protected)
        function y=apply_(this, x)
            % Reimplemented from parent class :class:`Cost`.
            if this.memoizeOpts.computeRes
                Res = this.memoize('computeRes', @this.computeRes_, x) ;
            else
                Res = this.computeRes_(x) ;
            end
            Res = abs(Res{1}) ;

            % method
            switch this.method
                case 'Andrews'
                    cond = Res < (pi*this.mu) ;
                    y = this.mu^2*( ...
                        cond.*(1-cos(Res)) ...
                        + (1-cond).*2) ;
                case 'Beaton-Tukey'
                    cond = Res < this.mu ;
                    y = this.mu^2/2*( ...
                        cond.*(1-(Res/this.mu).^2).^3 ...
                        + (1-cond)) ;
                case 'Cauchy'
                    y = this.mu^2/2*log(1+(Res/this.mu).^2) ;
                case 'Fair'
                    y = this.mu^2*(Res/this.mu-log(1+Res/this.mu)) ;
                case 'Huber'
                    cond = Res < this.mu ;
                    y = cond.*Res.^2/2 ...
                        + (1-cond).*(this.mu*Res-this.mu^2/2) ;
                case 'Logistic'
                    y = this.mu^2*log(cosh(Res/this.mu)) ;
                case 'Talwar-Hinich'
                    cond = Res < this.mu ;
                    y = cond.*Res.^2/2 ...
                        + (1-cond).*this.mu^2/2 ;
                case 'Welsch-Dennis'
                    y = this.mu^2/2*(1-exp(-(Res/this.mu).^2)) ;
            end
            y = sum(y(:)) ;
        end
        
        
        function g=applyGrad_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            if this.memoizeOpts.computeRes
                Res = this.memoize('computeRes', @this.computeRes_, x) ;
            else
                Res = this.computeRes_(x) ;
            end
            
            % Corresponding weights
            W = this.computeW_(x) ;
            
            % Composition
            if this.noise_model == 1 % Poisson noise
                Mx = Res{3} ;
                s = Res{2} ;
                Res = Res{1}.*W .* ...
                    -(this.eta.*(Mx+this.y)+2*this.var_0) ./ ...
                    (2.*s.*(this.eta.*Mx+this.var_0).^1.5) ;

            else % No noise model
                s = Res{2} ;
                Res = -Res{1}.*W./s ;
            end
            g = this.M.applyJacobianT_(Res, x) ;
        end
        
        
        % Function to compute and scale the residues
        function Res = computeRes_(this, x)
            Mx = this.M*x ;
            Res = this.y-Mx ;
            
            % Scaling by the factor s
            if ~isnumeric(this.flag_s)
                switch this.flag_s
                    case 'STD'
                        % standard deviation
                        s = 3 * std(Res(:)) ;
                    case 'MAD'
                        % median absolute deviation
                        s = 1.4826*median(abs(Res(:)-median(Res(:)))) ;
                end
            else
                s = this.flag_s ;
            end
            
            % Scaling according to the noise standard deviation if needed
            if this.noise_model == 1 % Poisson
                Res = Res./(s.*(this.eta.*Mx+this.var_0).^0.5) ;
                Res = {Res, s, Mx} ;
            else % No noise model
                Res = Res./s ;
                Res = {Res, s} ;
            end
        end
    end
    
    %% Methods (public)
    methods 
        % Function to compute weights corresponding to given residues
        function W = computeW_(this, x)
            if this.memoizeOpts.computeRes
                Res = this.memoize('computeRes', @this.computeRes_, x) ;
            else
                Res = this.computeRes_(x) ;
            end
            Res = abs(Res{1}) ;
            
            % method
            switch this.method
                case 'Andrews'
                    cond = Res < (pi*this.mu) ;
                    W = cond.*this.mu./Res.*sin(Res./this.mu) ;
                case 'Beaton-Tukey'
                    cond = Res < this.mu ;
                    W = cond.*(1-(Res/this.mu).^2).^2 ;
                case 'Cauchy'
                    W = 1./(1+(Res/this.mu).^2) ;
                case 'Fair'
                    W = 1./(1+abs(Res)/this.mu) ;
                case 'Huber'
                    cond = Res < this.mu ;
                    W = cond ...
                        + (1-cond).*(this.mu./abs(Res)) ;
                case 'Logistic'
                    W = this.mu./Res.*tanh(Res/this.mu) ;
                case 'Talwar-Hinich'
                    W = Res < this.mu ;
                case 'Welsch-Dennis'
                    W = exp(-(Res/this.mu).^2) ;
            end
        end
    end
end
