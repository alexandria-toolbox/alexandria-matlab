classdef MaximumLikelihoodVar < handle & VectorAutoRegression


    % Maximum likelihood vector autoregression, developed in section 11.1
    % 
    % Parameters:
    % -----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (4.11.1)
    % 
    % exogenous : matrix of size (n_obs,m), default = []
    %     exogenous variables, defined in (4.11.1)
    %
    % structural_identification : int, default = 2
    %     structural identification scheme, as defined in section 13.2
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % lags : int, default = 4
    %     number of lags, defined in (4.11.1)
    %
    % constant : bool, default = true
    %     if true, an intercept is included in the VAR model exogenous
    %
    % trend : bool, default = false
    %     if true, a linear trend is included in the VAR model exogenous
    %
    % quadratic_trend : bool, default = false
    %     if true, a quadratic trend is included in the VAR model exogenous
    % 
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    % 
    % verbose : bool, default = false
    %     if true, displays a progress bar  
    % 
    % 
    % Properties
    % ----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (4.11.1)
    % 
    % exogenous : matrix of size (n_obs,m)
    %     exogenous variables, defined in (4.11.1)
    %
    % structural_identification : int
    %     structural identification scheme, as defined in section 13.2
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % lags : int
    %     number of lags, defined in (4.11.1)
    %
    % constant : bool
    %     if true, an intercept is included in the VAR model exogenous
    %
    % trend : bool
    %     if true, a linear trend is included in the VAR model exogenous
    %
    % quadratic_trend : bool
    %     if true, a quadratic trend is included in the VAR model exogenous
    % 
    % credibility_level : float
    %     VAR model credibility level (between 0 and 1)
    % 
    % verbose : bool, default = false
    %     if true, displays a progress bar  
    % 
    % B : matrix of size (k,n)
    %     matrix of VAR coefficients, defined in (4.11.2)
    % 
    % Sigma : matrix of size (n,n)
    %     variance-covariance matrix of VAR residuals, defined in (4.11.1)
    % 
    % beta_estimates : matrix of size (k,n,3)
    %     estimates of VAR coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % Sigma_estimates : matrix of size (n,n)
    %     estimates of variance-covariance matrix of VAR residuals
    %
    % H :  matrix of size (n,n)
    %     structural identification matrix, defined in (4.13.5)
    %
    % Gamma : matrix of size (1,n)
    %     diagonal of structural shock variance matrix, defined in section 13.2
    %
    % Gamma_estimates : matrix of size (1,n)
    %     estimates of structural shock variance matrix, defined in section 13.2
    %
    % steady_state_estimates : matrix of size (T,n,3)
    %     estimates of steady-state, defined in (4.12.30)
    %
    % fitted_estimates : matrix of size (T,n,3)
    %     estimates of in-sample fit, defined in (4.11.2)
    %
    % residual_estimates : matrix of size (T,n,3)
    %     estimates of in-sample residuals, defined in (4.11.2)
    %
    % structural_shocks_estimates : matrix of size (T,n,3)
    %     estimates of in-sample structural shocks, defined in definition 13.1
    %
    % insample_evaluation : struct
    %     in-sample evaluation criteria, defined in (4.13.15)-(4.13.17)
    %
    % forecast_estimates : matrix of size (f_periods,n,3)
    %     forecast estimates, defined in (4.13.12) and (4.13.13)
    %     page 1: median, page 2: lower bound, page 3: upper bound
    %
    % forecast_evaluation_criteria : struct
    %     forecast evaluation criteria, defined in (4.13.18)-(4.13.19)
    %
    % irf_estimates : matrix of size (n,n,irf_periods,3)
    %     estimates of impulse response function, defined in (4.13.1) or (4.13.9)
    %     page 1: median, page 2: lower bound, page 3: upper bound    
    %
    % exo_irf_estimates : matrix of size (n,m,irf_periods,3)
    %     estimates of exogenous impulse response function, if any exogenous variable
    %     page 1: median, page 2: lower bound, page 3: upper bound
    %
    % fevd_estimates : matrix of size (n,n,fevd_periods,3)
    %     estimates of forecast error variance decomposition, defined in (4.13.31)
    %     page 1: median, page 2: lower bound, page 3: upper bound 
    %
    % hd_estimates : matrix of size (n,n,T,3)
    %     estimates of historical decomposition, defined in (4.13.35)
    %     page 1: median, page 2: lower bound, page 3: upper bound 
    % 
    % Y : matrix of size (T,n)
    %     matrix of in-sample endogenous variables, defined in (4.11.3)
    % 
    % Z : matrix of size (T,m)
    %     matrix of in-sample endogenous variables, defined in (4.11.3)
    % 
    % X : matrix of size (T,k)
    %     matrix of exogenous and lagged regressors, defined in (4.11.3)
    % 
    % n : int
    %     number of endogenous variables, defined in (4.11.1)
    % 
    % m : int
    %     number of exogenous variables, defined in (4.11.1)
    % 
    % p : int
    %     number of lags, defined in (4.11.1)
    % 
    % T : int
    %     number of sample periods, defined in (4.11.1)
    % 
    % k : int
    %     number of VAR coefficients in each equation, defined in (4.11.1)
    % 
    % q : int
    %     total number of VAR coefficients, defined in (4.11.1)
    %
    %
    % Methods
    % ----------
    % estimate
    % insample_fit
    % forecast
    % forecast_evaluation
    % impulse_response_function
    % forecast_error_variance_decomposition
    % historical_decomposition


    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess= protected)
        endogenous
        exogenous
        structural_identification
        lags
        constant
        trend
        quadratic_trend
        credibility_level
        verbose
        B
        Sigma
        beta_estimates
        Sigma_estimates
        H
        Gamma
        Gamma_estimates
        steady_state_estimates
        fitted_estimates
        residual_estimates
        structural_shock_estimates
        insample_evaluation
        forecast_estimates
        forecast_evaluation_criteria
        irf_estimates
        exo_irf_estimates
        fevd_estimates
        hd_estimates
    end    


    properties (GetAccess = private, SetAccess = private)
        inv_H
    end

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    
    
    methods (Access = public)
 
        
        function [self] = MaximumLikelihoodVar(endogenous, varargin)
            
            % constructor for the MaximumLikelihoodVar class
            
            % allow for optional arguments
            parser = inputParser;
            default_exogenous = [];
            default_structural_identification = 2;
            default_lags = 4;
            default_constant = true;
            default_trend = false;
            default_quadratic_trend = false; 
            default_credibility_level = 0.95;
            default_verbose = false;
            addRequired(parser,'endogenous');
            addParameter(parser,'exogenous', default_exogenous);
            addParameter(parser,'structural_identification', default_structural_identification);
            addParameter(parser,'lags', default_lags);
            addParameter(parser, 'constant', default_constant);
            addParameter(parser, 'trend', default_trend);
            addParameter(parser, 'quadratic_trend', default_quadratic_trend);            
            addParameter(parser,'credibility_level', default_credibility_level);
            addParameter(parser,'verbose', default_verbose);
            parse(parser, endogenous, varargin{:});
            self.endogenous = endogenous;
            self.exogenous = parser.Results.exogenous;
            self.structural_identification = parser.Results.structural_identification;
            self.lags = parser.Results.lags;
            self.constant = parser.Results.constant;
            self.trend = parser.Results.trend;
            self.quadratic_trend = parser.Results.quadratic_trend;            
            self.credibility_level = parser.Results.credibility_level;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
        end
        
        
        function estimate(self)
            
            % estimate()
            % estimates parameters B and Sigma of VAR model
            %
            % parameters:
            % none
            %
            % returns:
            % none

            % fit to obtain maximum likelihood estimates
            self.fit();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
            % estimate structural VAR
            self.make_structural_identification();
        end        

        
        function insample_fit(self)
            
            % fit_and_residuals()
            % estimates of in-sample fit elements (fit, residuals, in-sample criteria, steady-state, structural shocks)
            %
            % parameters:
            % none
            %
            % returns:
            % none
            
            % compute steady-state
            self.steady_state();
            % obtain fit and residuals
            self.fitted_and_residuals();
            % obtain in-sample evaluation criteria
            self.insample_criteria();     
        end
        
        
        function [forecast_estimates] = forecast(self, h, credibility_level, varargin)
            
            % forecast(self, h, credibility_level, varargin)
            % predictions for the maximum likelihood VAR model using (4.13.12)
            %
            % parameters:
            % h : int
            %     number of forecast periods
            % credibility_level : float
            %     credibility level for predictions (between 0 and 1)          
            % varargin : either zero or one additional argument, matrix of dimension (h, n_exo)
            %     no argument unless the model includes exogenous other than constant, trend and quadratic trend
            %     if argument, n_exo is the number of additional exogenous variables          
            %
            % returns:
            % forecast_estimates : matrix of size (periods, n, 3)
            %    posterior estimates for predictions
            %    page 1: median; page 2: interval lower bound;
            %    page 3: interval upper bound         

            % optional argument
            if isempty(varargin)
                Z_p = [];
            else
                Z_p = varargin{1};
            end  
            % unpack
            Y = self.Y;
            B = self.B;
            Sigma = self.Sigma;
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            n = self.n;
            p = self.p;
            T = self.T;
            exogenous = self.exogenous;
            % make regressors
            [Z_p, Y] = vu.make_forecast_regressors(Z_p, Y, h, p, T, exogenous, constant, trend, quadratic_trend);
            % obtain point estimates
            [forecasts] = vu.linear_forecast(B, h, Z_p, Y, n);
            % obtain confidence intervals
            [lower_bound, upper_bound] = self.maximum_likelihood_forecast_credibility( ...
            forecasts, credibility_level, B, Sigma, h, n, p);
            % concatenate
            forecast_estimates = cat(3,forecasts,lower_bound,upper_bound);
            self.forecast_estimates = forecast_estimates;
        end
        
        
        function forecast_evaluation(self, Y)
            
            % forecast_evaluation(Y)
            % forecast evaluation criteria for the maximum likelihood VAR model
            % 
            % parameters:
            % Y : ndarray of shape (perods,n)
            %     array of realised values for forecast evaluation
            % 
            % returns:
            %     none

            % recover forecasts
            Y_hat = self.forecast_estimates(:,:,1);
            % obtain forecast evaluation criteria from equations (4.13.18) and (4.13.19)
            forecast_evaluation_criteria = vu.forecast_evaluation_criteria(Y_hat, Y);
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;
        end


        function [irf_estimates exo_irf_estimates] = impulse_response_function(self, h, credibility_level)
            
            % impulse_response_function(h, credibility_level)
            % impulse response functions, as defined in (4.13.1)-(4.13.9)
            %
            % parameters:
            % h : int
            %     number of IRF periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            %    
            % returns:
            % irf_estimates : ndarray of shape (n,n,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound
            % exo_irf_estimates : ndarray of shape (n,m,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound
            
            % get regular impulse response funtion
            [irf mcmc_irf irf_exo mcmc_irf_exo] = self.make_impulse_response_function(h);    
            % get structural impulse response function
            [structural_irf mcmc_structural_irf] = self.make_structural_impulse_response_function(irf, mcmc_irf, h);
            % obtain posterior estimates
            [irf_estimates exo_irf_estimates] = self.irf_posterior_estimates(credibility_level, irf, mcmc_irf, ...
                                                irf_exo, mcmc_irf_exo, structural_irf, mcmc_structural_irf);
            self.irf_estimates = irf_estimates;
            self.exo_irf_estimates = exo_irf_estimates;           
        end


        function [fevd_estimates] = forecast_error_variance_decomposition(self, h, credibility_level)
            
            % forecast_error_variance_decomposition(self, h, credibility_level)
            % forecast error variance decomposition, as defined in (4.13.31)
            % 
            % parameters:
            % h : int
            %     number of FEVD periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands            
            %   
            % returns:
            % fevd_estimates : matrix of size (n,n,h)
            %     dimensions are variable, shock, period
    
            % get forecast error variance decomposition
            [fevd mcmc_fevd] = self.make_forecast_error_variance_decomposition(h);
            % obtain posterior estimates
            fevd_estimates = self.fevd_posterior_estimates(credibility_level, fevd, mcmc_fevd);
            self.fevd_estimates = fevd_estimates;
            if self.verbose
                cu.progress_bar_complete('Forecast error variance decomposition:');
            end
        end


        function [hd_estimates] = historical_decomposition(self, credibility_level)

            % [hd_estimates] = historical_decomposition(self, credibility_level)
            % historical decomposition, as defined in (4.13.34)-(4.13-36)
            % 
            % parameters:
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % hd_estimates : matrix of size (n,n,T,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound

            % get historical decomposition
            [hd mcmc_hd] = self.make_historical_decomposition();
            % obtain posterior estimates
            hd_estimates = self.hd_posterior_estimates(credibility_level, hd, mcmc_hd);
            self.hd_estimates = hd_estimates;
            if self.verbose
                cu.progress_bar_complete('Historical decomposition:');
            end
        end  


    end   
    
        
    methods (Access = protected, Hidden = true)        
        

        function fit(self)
            
            % estimates B_hat and Sigma_hat from (4.11.9)

            % unpack        
            Y = self.Y;
            X = self.X;
            XX = self.XX;
            XY = self.XY;
            T = self.T;
            verbose = self.verbose;
            % estimate B_hat and Sigma_hat
            [B, Sigma] = self.ols_var(Y, X, XX, XY, T);
            if verbose
                cu.progress_bar_complete('VAR parameters:');
            end
            self.B = B;
            self.Sigma = Sigma;
        end
        
        
        function parameter_estimates(self)
            
            % estimates and intervals from Normal distribution in (4.11.10)

            % unpack
            B = self.B;
            Sigma = self.Sigma;
            credibility_level = self.credibility_level;
            XX = self.XX;
            k = self.k;
            n = self.n;
            % get coefficients variance Q
            Q = kron(Sigma, la.invert_spd_matrix(XX));
            Qii = reshape(sqrt(diag(Q)),[k n]);
            % critical value of Normal distribution for credibility level
            Z = su.normal_icdf((1 + credibility_level) / 2);
            % initiate storage: 4 dimensions: median, standard deviation, lower bound, upper bound
            beta_estimates = zeros(k,n,4); 
            beta_estimates(:,:,1) = B;
            beta_estimates(:,:,2) = B - Z * Qii;
            beta_estimates(:,:,3) = B + Z * Qii; 
            beta_estimates(:,:,4) = Qii;
            % save as attributes
            self.beta_estimates = beta_estimates;
            self.Sigma_estimates = self.Sigma;
        end
        
        
        function make_structural_identification(self)
            
            % computes structural VAR estimates
            
            % unpack
            n = self.n;
            structural_identification = self.structural_identification;
            Sigma = self.Sigma;
            verbose = self.verbose;
            % estimate Cholesky SVAR, if selected
            if structural_identification == 2
                H = la.cholesky_nspd(Sigma);
                inv_H = la.invert_lower_triangular_matrix(H);
                self.H = H;
                self.inv_H = inv_H;
                self.Gamma = ones(1,n); 
                self.Gamma_estimates = ones(1,n); 
                if verbose
                    cu.progress_bar_complete('SVAR parameters:');
                end
            % estimate triangular SVAR, if selected
            elseif structural_identification == 3
                [H, Gamma] = la.triangular_factorization(Sigma);
                inv_H = la.invert_lower_triangular_matrix(H);
                self.H = H;
                self.inv_H = inv_H;
                self.Gamma = Gamma; 
                self.Gamma_estimates = Gamma;
                if verbose
                    cu.progress_bar_complete('SVAR parameters:');
                end
            end
        end

        
        function steady_state(self)   

            % computes steady-state for the VAR model
            
            % point estimates
            Z = self.Z;
            B = self.B;
            Sigma = self.Sigma;
            n = self.n;
            m = self.m;
            p = self.p;
            T = self.T;
            k = self.k;
            q = self.q;
            XX = self.XX;
            ss = vu.steady_state(Z, B, n, m, p, T);
            % simulated coefficients for credibility interval
            mcmc_B = vu.ols_var_mcmc_beta(B, Sigma, XX, k, n, q);
            mcmc_ss = zeros(T,n,500);
            for j=1:500
                mcmc_ss(:,:,j) = vu.steady_state(Z, mcmc_B(:,:,j), n, m, p, T);
            end
            if self.verbose 
                cu.progress_bar_complete('Steady-state:')
            end
            ss_estimates = vu.posterior_estimates(mcmc_ss, self.credibility_level);
            ss_estimates(:,:,1) = ss;
            self.steady_state_estimates = ss_estimates;
        end


        function fitted_and_residuals(self)
            
            % computes fit and residuals

            % point estimates
            Y = self.Y;            
            X = self.X;            
            B = self.B;  
            Sigma = self.Sigma;
            n = self.n;
            T = self.T;
            k = self.k;
            q = self.q;
            XX = self.XX;
            structural_identification = self.structural_identification;
            % simulated coefficients for credibility interval
            mcmc_B = vu.ols_var_mcmc_beta(B, Sigma, XX, k, n, q);
            mcmc_fitted = zeros(T,n,500);
            mcmc_residuals = zeros(T,n,500);
            for j=1:500
                [E, Y_hat] = vu.fit_and_residuals(Y, X, mcmc_B(:,:,j)); 
                mcmc_fitted(:,:,j) = Y_hat;
                mcmc_residuals(:,:,j) = E;
            end
            if self.verbose 
                cu.progress_bar_complete('Fitted and residual:')
            end
            fitted_estimates = vu.posterior_estimates(mcmc_fitted, self.credibility_level);
            residual_estimates = vu.posterior_estimates(mcmc_residuals, self.credibility_level);
            [E, Y_hat] = vu.fit_and_residuals(Y, X, B);  
            fitted_estimates(:,:,1) = Y_hat;
            residual_estimates(:,:,1) = E;
            if ismember(structural_identification, [2,3])
                inv_H = self.inv_H;
                mcmc_structural_shocks = zeros(T,n,500);
                for j=1:500
                    Xi = vu.structural_shocks(mcmc_residuals(:,:,j), inv_H);
                    mcmc_structural_shocks(:,:,j) = Xi;
                end
                if self.verbose
                    cu.progress_bar_complete('Structural shocks:')
                end
                structural_shock_estimates = vu.posterior_estimates(mcmc_structural_shocks, self.credibility_level);
                [E, Y_hat] = vu.fit_and_residuals(Y, X, B);
                Xi = vu.structural_shocks(E, inv_H);
                structural_shock_estimates(:,:,1) = Xi;
            else
                structural_shock_estimates = [];
            end
            self.fitted_estimates = fitted_estimates;
            self.residual_estimates = residual_estimates;
            self.structural_shock_estimates = structural_shock_estimates;
        end

   
        function insample_criteria(self)
            
            % computes in-sample evluation criteria

            Y = self.Y;   
            E = self.residual_estimates(:,:,1); 
            T = self.T;   
            k = self.k;          
            Sigma = self.Sigma;            
            q = self.q;       
            % estimate general criteria (SSR, R2, adj-R2)
            insample_evaluation_1 = vu.insample_evaluation_criteria(Y, E, T, k);       
            % estimate criteria specific to maximum likelihood VAR (aic, bic, hq)
            insample_evaluation_2 = self.maximum_likelihood_evaluation_criteria(Sigma, q, T);          
            % merge in-sample criteria
            insample_evaluation = iu.concatenate_structures(insample_evaluation_1, insample_evaluation_2);          
            self.insample_evaluation = insample_evaluation;         
        end

        
        function [insample_evaluation] = maximum_likelihood_evaluation_criteria(self, Sigma, q, T)      

            % computes aic, bic and hq for maximum likelihood VAR
            
            log_det_Sigma = log(det(Sigma));
            aic = 2 * q / T + log_det_Sigma;
            bic = q * log(T) / T + log_det_Sigma;
            hq = 2 * q * log(log(T)) / T + log_det_Sigma;
            insample_evaluation = struct;
            insample_evaluation.aic = aic;
            insample_evaluation.bic = bic;
            insample_evaluation.hq = hq;     
        end
    
        
        function [lower_bound, upper_bound] = maximum_likelihood_forecast_credibility(self, forecasts, credibility_level, B, Sigma, periods, n, p)
    
            % create forecast credibilty intervals for the maximum likelihood VAR
    
            lower_bound = zeros(periods,n);
            upper_bound = zeros(periods,n);
            Q_h = zeros(n,n);
            Z = su.normal_icdf((1 + credibility_level) / 2);
            irf = vu.impulse_response_function(B, n, p, periods);
            for h=1:periods
                y_h = forecasts(h,:);
                Phi_h = irf(:,:,h);
                Q_h = Q_h + Phi_h * Sigma * Phi_h';
                s_h = sqrt(diag(Q_h))';
                lower_bound(h,:) = y_h - Z * s_h;
                upper_bound(h,:) = y_h + Z * s_h;
            end          
        end        


        function [irf mcmc_irf irf_exo mcmc_irf_exo] = make_impulse_response_function(self, h)
            
            % point estimates and simulations for impulse response functions
            
            % simulated coefficients for credibility interval
            mcmc_B = vu.ols_var_mcmc_beta(self.B, self.Sigma, self.XX, self.k, self.n, self.q);
            % impulse response functions, point estimates
            irf = vu.impulse_response_function(self.B, self.n, self.p, h);
            % impulse response functions, simulated values
            mcmc_irf = zeros(self.n, self.n, h, 500);
            for i=1:500
                mcmc_irf(:,:,:,i) = vu.impulse_response_function(mcmc_B(:,:,i), self.n, self.p, h);
            end
            if self.verbose 
                cu.progress_bar_complete('Impulse response function:')
            end
            % exogenous impulse response functions, point estimates
            if ~isempty(self.exogenous)
                r = size(self.exogenous,2);
                irf_exo = vu.exogenous_impulse_response_function(self.B, self.n, self.m, r, self.p, h);
                % impulse response functions, simulated values
                mcmc_irf_exo = zeros(self.n, r, h, 500);
                for i=1:500
                    mcmc_irf_exo(:,:,:,i) = vu.exogenous_impulse_response_function(mcmc_B(:,:,i), self.n, self.m, r, self.p, h);
                end
                if self.verbose  
                    cu.progress_bar_complete('Exogenous impulse response function:');
                end
            else
                irf_exo = [];
                mcmc_irf_exo = [];
            end
        end


        function [structural_irf mcmc_structural_irf] = make_structural_impulse_response_function(self, irf, mcmc_irf, h)
            
            % structural impulse response function
    
            % get structural impulse response function
            if self.structural_identification == 1
                structural_irf = [];
                mcmc_structural_irf = [];
            else
                structural_irf = vu.structural_impulse_response_function(irf, self.H, self.n);
                mcmc_structural_irf = zeros(self.n, self.n, h, 500);
                for i=1:500
                    mcmc_structural_irf(:,:,:,i) = vu.structural_impulse_response_function(mcmc_irf(:,:,:,i), self.H, self.n);
                end
                if self.verbose
                    cu.progress_bar_complete('Structural impulse response function:');
                end
            end
        end


        function [irf_estimates exo_irf_estimates] = irf_posterior_estimates(self, credibility_level, irf, ...
                                                     mcmc_irf, irf_exo, mcmc_irf_exo, structural_irf, mcmc_structural_irf)
            
            % posterior estimates of impulse response function
            
            % posterior estimates for endogenous
            if self.structural_identification == 1
                irf_estimates = vu.posterior_estimates_3d(mcmc_irf, credibility_level);
                irf_estimates(:,:,:,1) = irf;
            else
                irf_estimates = vu.posterior_estimates_3d(mcmc_structural_irf, credibility_level);
                irf_estimates(:,:,:,1) = structural_irf;
            end
            % posterior estimates for exogenous
            if ~isempty(self.exogenous)
                exo_irf_estimates = vu.posterior_estimates_3d(mcmc_irf_exo, credibility_level);
            else
                exo_irf_estimates = [];
            end
        end
        

        function [fevd mcmc_fevd] = make_forecast_error_variance_decomposition(self, h)
            
            % forecast error variance decomposition
            
            % if no structural identification, fevd is not computed
            if self.structural_identification == 1
                fevd = [];
                mcmc_fevd = [];
            else
                % fevd, point estimates
                irf = vu.impulse_response_function(self.B, self.n, self.p, h);
                structural_irf = vu.structural_impulse_response_function(irf, self.H, self.n);
                fevd = vu.forecast_error_variance_decomposition(structural_irf, self.Gamma, self.n, h);
                % fevd, simulated values
                mcmc_B = vu.ols_var_mcmc_beta(self.B, self.Sigma, self.XX, self.k, self.n, self.q);
                mcmc_fevd = zeros(self.n, self.n, h, 500);
                for i=1:500
                    irf = vu.impulse_response_function(mcmc_B(:,:,i), self.n, self.p, h);
                    structural_irf = vu.structural_impulse_response_function(irf, self.H, self.n); 
                    mcmc_fevd(:,:,:,i) = vu.forecast_error_variance_decomposition(structural_irf, self.Gamma, self.n, h);
                end
            end
        end


        function [fevd_estimates] = fevd_posterior_estimates(self, credibility_level, fevd, mcmc_fevd)
    
            % posterior estimates of forecast error variance decomposition
            
            % posterior estimates
            if self.structural_identification == 1
                fevd_estimates = [];
            else
                fevd_estimates = vu.posterior_estimates_3d(mcmc_fevd, credibility_level);
                fevd_estimates(:,:,:,1) = fevd;
            end
        end


        function [hd mcmc_hd] = make_historical_decomposition(self)
            
            % historical decomposition
            
            % if no structural identification, historical decomposition is not computed
            if self.structural_identification == 1
                hd = [];
                mcmc_hd = [];
            else
                % hd, point estimates
                irf = vu.impulse_response_function(self.B, self.n, self.p, self.T);
                structural_irf = vu.structural_impulse_response_function(irf, self.H, self.n);
                [E ~] = vu.fit_and_residuals(self.Y, self.X, self.B);
                structural_shocks = vu.structural_shocks(E, self.inv_H);    
                hd = vu.historical_decomposition(structural_irf, structural_shocks, self.n, self.T); 
                % hd, simulated values
                mcmc_B = vu.ols_var_mcmc_beta(self.B, self.Sigma, self.XX, self.k, self.n, self.q);
                mcmc_hd = zeros(self.n, self.n, self.T, 500);
                for i=1:500
                    irf = vu.impulse_response_function(mcmc_B(:,:,i), self.n, self.p, self.T);
                    structural_irf = vu.structural_impulse_response_function(irf, self.H, self.n);
                    [E ~] = vu.fit_and_residuals(self.Y, self.X, mcmc_B(:,:,i));
                    structural_shocks = vu.structural_shocks(E, self.inv_H);           
                    mcmc_hd(:,:,:,i) = vu.historical_decomposition(structural_irf, structural_shocks, self.n, self.T);
                end
            end
        end


        function [hd_estimates] = hd_posterior_estimates(self, credibility_level, hd, mcmc_hd)
    
            % posterior estimates of historical decomposition
            
            % posterior estimates
            if self.structural_identification == 1
                hd_estimates = [];
            else
                hd_estimates = vu.posterior_estimates_3d(mcmc_hd, credibility_level);
                hd_estimates(:,:,:,1) = hd;
            end
        end

        
    end
    
    
end

