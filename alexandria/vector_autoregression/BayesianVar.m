classdef BayesianVar < handle
    
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        delta
        s
        Y_sum
        X_sum
        Y_obs
        X_obs
        Y_lrp
        X_lrp
        Y_d
        X_d
        T_d
        steady_state_estimates      
        fitted_estimates
        residual_estimates
        structural_shock_estimates
        insample_evaluation
        mcmc_structural_shocks
        mcmc_forecast
        forecast_estimates
        forecast_evaluation_criteria
        mcmc_irf
        mcmc_irf_exo
        mcmc_structural_irf
        irf_estimates
        exo_irf_estimates
        mcmc_fevd
        fevd_estimates
        mcmc_hd
        hd_estimates
        mcmc_conditional_forecast
        conditional_forecast_estimates
        H_estimates
        Gamma_estimates
    end    
    
    properties (GetAccess = protected, SetAccess = protected)
        XX_d
        XY_d
        YY_d
        mcmc_inv_H
        svar_index
        secondary_constrained_coefficients_table
    end 

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public) 
    
    
        function self = BayesianVar()
        end
        
        
        function insample_fit(self)
            
            % insample_fit(self)
            % generates in-sample fit and residuals along with evaluation criteria
            %
            % parameters:
            % none
            %
            % returns:
            % none
            
            % compute steady-state
            self.steady_state();
            % compute fitted and residuals
            self.fitted_and_residual();
            % compute in-sample criteria
            self.insample_criteria();
        end            

        
        function [forecast_estimates] = forecast(self, h, credibility_level, varargin)
            
            % [forecast_estimates] = forecast(self, h, credibility_level, varargin)
            % estimates forecasts for the Bayesian VAR model, using algorithm 13.4
            % 
            % parameters:
            % h : int
            %     number of forecast periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % varargin : either zero or one additional argument, matrix of dimension (h, n_exo)
            %     no argument unless the model includes exogenous other than constant, trend and quadratic trend
            %     if argument, n_exo is the number of additional exogenous variables
            % 
            % returns:
            % forecast_estimates : matrix of dimension (h,n,3)
            %     page 1: median; page 2: interval lower bound; page 3: interval upper bound            
            
            % optional argument
            if isempty(varargin)
                Z_p = [];
            else
                Z_p = varargin{1};
            end  
            % get forecast
            self.make_forecast(h, Z_p);
            % obtain posterior estimates
            self.forecast_posterior_estimates(credibility_level);
            forecast_estimates = self.forecast_estimates;        
        end      
        

        function [forecast_evaluation_criteria] = forecast_evaluation(self, Y)
            
            % [forecast_evaluation_criteria] = forecast_evaluation(self, Y)
            % forecast evaluation criteria for the Bayesian VAR model, as defined in (4.13.18)-(4.13.22)
            % 
            % parameters:
            % Y : matrix of size(h,n)
            %     matrix of realised values for forecast evaluation, h being the number of forecast periods
            % 
            % returns:
            % forecast_evaluation_criteria : structure
            %     structure with criteria name as keys and corresponding number as value       
            
            % unpack
            Y_hat = self.forecast_estimates(:,:,1);
            mcmc_forecast = self.mcmc_forecast;
            % obtain regular forecast evaluation criteria from equations (4.13.18) and (4.13.19)
            standard_evaluation_criteria = vu.forecast_evaluation_criteria(Y_hat, Y);
            % obtain Bayesian forecast evaluation criteria from equations (4.13.21) and (4.13.22)
            bayesian_evaluation_criteria = vu.bayesian_forecast_evaluation_criteria(mcmc_forecast, Y);
            % merge structures
            forecast_evaluation_criteria = iu.concatenate_structures(standard_evaluation_criteria, bayesian_evaluation_criteria);
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;     
        end    
        
        
        function [irf_estimates exo_irf_estimates] = impulse_response_function(self, h, credibility_level)

            % [irf_estimates exo_irf_estimates] = impulse_response_function(h, credibility_level)
            % impulse response functions, as defined in (4.13.1)-(4.13.9)
            % 
            % parameters:
            % h : int
            %     number of IRF periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % irf_estimates : matrix of size (n,n,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound
            % exo_irf_estimates : matrix of size (n,n,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound

            % get regular impulse response funtion
            self.make_impulse_response_function(h);
            % get exogenous impuse response function
            self.make_exogenous_impulse_response_function(h);
            % get structural impulse response function
            self.make_structural_impulse_response_function(h);
            % obtain posterior estimates
            self.irf_posterior_estimates(credibility_level);
            irf_estimates = self.irf_estimates;
            exo_irf_estimates = self.exo_irf_estimates;
        end    
        
        
        function [fevd_estimates] = forecast_error_variance_decomposition(self, h, credibility_level)

            % [fevd_estimates] = forecast_error_variance_decomposition(self, h, credibility_level)
            % forecast error variance decomposition, as defined in (4.13.31)
            % 
            % parameters:
            % h : int
            %     number of FEVD periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % fevd_estimates : matrix of size (n,n,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound

            % get forecast error variance decomposition
            self.make_forecast_error_variance_decomposition(h);
            % obtain posterior estimates
            self.fevd_posterior_estimates(credibility_level);
            fevd_estimates = self.fevd_estimates;
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
            self.make_historical_decomposition();
            % obtain posterior estimates
            self.hd_posterior_estimates(credibility_level);
            hd_estimates = self.hd_estimates;
        end    
        

        function [conditional_forecast_estimates] = conditional_forecast(self, h, credibility_level, conditions, shocks, conditional_forecast_type, varargin)

            % conditional_forecast(h, credibility_level, conditions, shocks, conditional_forecast_type, varargin)
            % estimates conditional forecasts for the Bayesian VAR model, using algorithms 14.1 and 14.2
            % 
            % parameters:
            % h : int
            %     number of forecast periods           
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % conditions: matrix of size (n_conditions,4)
            %     table defining conditions (column 1: variable, column 2: period, column 3: mean, column 4: variance)   
            % shocks: empty matrix or matrix of size (n,1)
            %     vector defining shocks generating the conditions; should be empty if conditional_forecast_type = 1  
            % conditional_forecast_type : int
            %     conditional forecast type (1 = agnostic, 2 = structural)
            % varargin : either zero or one additional argument, matrix of dimension (h, n_exo)
            %     no argument unless the model includes exogenous other than constant, trend and quadratic trend
            %     if argument, n_exo is the number of additional exogenous variables            
            % 
            % returns:
            % conditional_forecast_estimates : matrix of size (h,n,3)
            %     page 1: median; page 2: interval lower bound; page 3: interval upper bound

            % optional argument
            if isempty(varargin)
                Z_p = [];
            else
                Z_p = varargin{1};
            end 
            % if conditional forecast type is agnostic
            if conditional_forecast_type == 1
                % get conditional forecasts
                self.make_conditional_forecast(h, conditions, Z_p);
            % if instead conditional forecast type is structural
            elseif conditional_forecast_type == 2
                % establish type of shocks
                shock_type = self.check_shock_type(h, conditions, shocks);
                % get structural conditional forecasts
                self.make_structural_conditional_forecast(h, conditions, shocks, Z_p, shock_type);
            end
            % obtain posterior estimates
            self.conditional_forecast_posterior_estimates(credibility_level);
            conditional_forecast_estimates = self.conditional_forecast_estimates;
        end  
  
        
    end
    
    
    methods (Access = protected, Hidden = true) 
        

        function make_delta(self)
            
            if isscalar(self.ar_coefficients)
                ar_coefficients = repmat(self.ar_coefficients, [self.n 1]);
            else
                ar_coefficients = self.ar_coefficients;
            end
            self.delta = ar_coefficients;
        end
        
        
        function individual_ar_variances(self)
            
            s = zeros(self.n,1);
            for i=1:self.n
                ar = MaximumLikelihoodVar(self.endogenous(:,i), 'lags', self.lags);
                ar.estimate();
                s(i) = ar.Sigma;
            end
            self.s = s;
        end
        
        
        function dummy_extensions(self)
            
            % sums-of-coefficients
            [Y_sum X_sum] = vu.sums_of_coefficients_extension(self.sums_of_coefficients, ...
                                                              self.pi5, self.Y, self.n, self.m, self.p);
            % dummy initial observation
            [Y_obs X_obs] = vu.dummy_initial_observation_extension(self.dummy_initial_observation, ...
                                                                   self.pi6, self.Y, self.X, self.n, self.m, self.p);
            % long run prior
            [Y_lrp X_lrp] = vu.long_run_prior_extension(self.long_run_prior, self.pi7, self.J, ...
                                                        self.Y, self.n, self.m, self.p);
            % concatenate to Y, X, update T as in 4.11.62
            Y_d = [Y_sum;Y_obs;Y_lrp;self.Y];
            X_d = [X_sum;X_obs;X_lrp;self.X];
            T_d = size([Y_sum;Y_obs;Y_lrp],1) + self.T;
            self.Y_sum = Y_sum;
            self.X_sum = X_sum;
            self.Y_obs = Y_obs;
            self.X_obs = X_obs;
            self.Y_lrp = Y_lrp;
            self.X_lrp = X_lrp;
            self.XX_d = X_d' * X_d;
            self.XY_d = X_d' * Y_d;
            self.YY_d = Y_d' * Y_d;
            self.Y_d = Y_d;
            self.X_d = X_d;
            self.T_d = T_d;
        end
        
        
        function make_constrained_coefficients(self)
            
            % apply only if constrained coefficients is activated
            if self.constrained_coefficients
                % unpack
                n = self.n;
                m = self.m;
                k = self.k;
                lags = self.lags;
                constant = self.constant;
                trend = self.trend;
                quadratic_trend = self.quadratic_trend;
                constrained_coefficients_table = self.constrained_coefficients_table;
                secondary_constrained_coefficients_table = vu.rework_constraint_table(constrained_coefficients_table, lags);
                self.secondary_constrained_coefficients_table = secondary_constrained_coefficients_table;
                % reshape b and V as k * n arrays
                B = reshape(self.b, [k n]);
                V = reshape(self.V, [k n]);
                % get updated parameters with constrained coefficients applied
                [new_b new_V] = vu.make_constrained_coefficients(B, V, n, m, k, lags, constant,...
                                trend, quadratic_trend, secondary_constrained_coefficients_table);
                self.b = new_b;
                self.V = new_V;
            end
        end
        

        function make_structural_identification(self)
            
            % SVAR by Choleski factorization, if selected
            if self.structural_identification == 2
                self.svar_by_choleski_factorization();
            % SVAR by triangular factorization, if selected
            elseif self.structural_identification == 3
                self.svar_by_triangular_factorization();
            % SVAR by restrictions, if selected
            elseif self.structural_identification == 4
                self.svar_by_restrictions();
            end
            % get posterior estimates
            if self.structural_identification ~= 1
                self.svar_estimates();  
            end
        end
        
        
        function svar_by_choleski_factorization(self)
            
            % compute H, Gamma, and inverse of H
            self.mcmc_H = self.mcmc_chol_Sigma;
            self.mcmc_Gamma = ones(self.iterations,self.n);
            % check if Sigma is attribute, which means prior is Minnesota
            if isprop(self,'Sigma')
                self.mcmc_inv_H = repmat(la.invert_lower_triangular_matrix(...
                                  self.mcmc_H(:,:,1)), [1 1 self.iterations]);
                if self.verbose
                    cu.progress_bar_complete('Structural identification:');  
                end
            % for any other prior, use mcmc draws
            else
                mcmc_inv_H = zeros(self.n,self.n,self.iterations);
                for j=1:self.iterations
                    mcmc_inv_H(:,:,j) = la.invert_lower_triangular_matrix(self.mcmc_H(:,:,j));
                    if self.verbose
                        cu.progress_bar(j, self.iterations, 'Structural identification:');
                    end
                end
                self.mcmc_inv_H = mcmc_inv_H;
            end 
            % create index correspondance for later applications
            self.svar_index = (1:self.iterations)';
        end
        
        
        function svar_by_triangular_factorization(self)
            
            % check if Sigma is attribute, which means prior is Minnesota
            if isprop(self,'Sigma')
                [H Gamma] = la.triangular_factorization(self.Sigma);
                mcmc_H = repmat(H,[1 1 self.iterations]);
                mcmc_inv_H = repmat(la.invert_lower_triangular_matrix(H),[1 1 self.iterations]);
                mcmc_Gamma = repmat(Gamma,[self.iterations 1]);
                if self.verbose
                    cu.progress_bar_complete('Structural identification:');
                end
            % for any other prior, use mcmc draws
            else
                mcmc_H = zeros(self.n,self.n,self.iterations);
                mcmc_inv_H = zeros(self.n,self.n,self.iterations);
                mcmc_Gamma = zeros(self.iterations,self.n);
                for j=1:self.iterations
                    [H Gamma] = la.triangular_factorization(self.mcmc_chol_Sigma(:,:,j), true);
                    mcmc_H(:,:,j) = H;
                    mcmc_inv_H(:,:,j) = la.invert_lower_triangular_matrix(H);
                    mcmc_Gamma(j,:) = Gamma;
                end
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Structural identification:');
                end
            end
            self.mcmc_H = mcmc_H;
            self.mcmc_inv_H = mcmc_inv_H;
            self.mcmc_Gamma = mcmc_Gamma;
            % create index correspondance for later applications
            self.svar_index = (1:self.iterations)';
        end

        
        function svar_by_restrictions(self)
            
            % initiate MCMC elements
            svar_index = zeros(self.iterations,1);
            mcmc_H = zeros(self.n,self.n,self.iterations);
            mcmc_inv_H = zeros(self.n,self.n,self.iterations);
            % create matrices of restriction checks
            [restriction_matrices max_irf_period] = vu.make_restriction_matrices(self.restriction_table, self.p);
            % make preliminary orthogonalised impulse response functions, if relevant
            [mcmc_irf] = vu.make_restriction_irf(self.mcmc_beta, self.mcmc_chol_Sigma, ...
                                                 self.iterations, self.n, self.p, max_irf_period);  
            % make preliminary structural shocks, if relevant
            [mcmc_shocks] = vu.make_restriction_shocks(self.mcmc_beta, self.mcmc_chol_Sigma, self.Y, self.X, ...
                                                       self.T, self.n, self.iterations, restriction_matrices);  
            % loop over iterations, until desired number of total iterations is obtained
            i = 1;
            while i <= self.iterations
                % select a random index in number of iterations
                j = rn.discrete_uniform(1, self.iterations);
                % make a random rotation matrix Q: if no zero restrictions, draw from uniform orthogonal distribution
                if isempty(restriction_matrices{1,1})
                    Q = rn.uniform_orthogonal(self.n);
                % if there are zero restrictions, use the zero uniform orthogonal distribution
                else
                    Q = rn.zero_uniform_orthogonal(self.n, restriction_matrices{1,1}, mcmc_irf(:,:,:,j));
                end
                % check restrictions: IRF, sign
                irf_sign_index = restriction_matrices{2,1};
                if ~isempty(irf_sign_index)
                    irf_sign_coefficients = restriction_matrices{2,2};
                    restriction_satisfied = vu.check_irf_sign(irf_sign_index, irf_sign_coefficients, mcmc_irf(:,:,:,j), Q);
                    if ~restriction_satisfied
                        continue
                    end
                end
                % check restrictions: IRF, magnitude
                irf_magnitude_index = restriction_matrices{3,1};
                if ~isempty(irf_magnitude_index)
                    irf_magnitude_coefficients = restriction_matrices{3,2};
                    restriction_satisfied = vu.check_irf_magnitude(irf_magnitude_index, irf_magnitude_coefficients, mcmc_irf(:,:,:,j), Q);
                    if ~restriction_satisfied
                        continue
                    end
                end
                % check restrictions: structural shocks, sign
                shock_sign_index = restriction_matrices{4,1};
                if ~isempty(shock_sign_index)
                    shock_sign_coefficients = restriction_matrices{4,2};
                    restriction_satisfied = vu.check_shock_sign(shock_sign_index, shock_sign_coefficients, mcmc_shocks(:,:,j), Q);
                    if ~restriction_satisfied
                        continue
                    end
                end
                % check restrictions: structural shocks, magnitude
                shock_magnitude_index = restriction_matrices{5,1};
                if ~isempty(shock_magnitude_index)
                    shock_magnitude_coefficients = restriction_matrices{5,2};
                    restriction_satisfied = vu.check_shock_magnitude(shock_magnitude_index, shock_magnitude_coefficients, mcmc_shocks(:,:,j), Q);
                    if ~restriction_satisfied
                        continue                 
                    end
                end
                % historical decomposition values if any of sign or magnitude restrictions apply
                history_sign_index = restriction_matrices{6,1};
                history_magnitude_index = restriction_matrices{7,1};
                if ~isempty(history_sign_index) || ~isempty(history_magnitude_index)
                    [irf shocks] = vu.make_restriction_irf_and_shocks(mcmc_irf(:,:,:,j), mcmc_shocks(:,:,j), Q, self.n);             
                end
                % check restrictions: historical decomposition, sign
                if ~isempty(history_sign_index)
                    history_sign_coefficients = restriction_matrices{6,2};
                    restriction_satisfied = vu.check_history_sign(history_sign_index, history_sign_coefficients, irf, shocks);
                    if ~restriction_satisfied 
                        continue                
                    end
                end
                % check restrictions: historical decomposition, magnitude
                if ~isempty(history_magnitude_index)
                    history_magnitude_coefficients = restriction_matrices{7,2};
                    restriction_satisfied = vu.check_history_magnitude(history_magnitude_index, history_magnitude_coefficients, irf, shocks);
                    if ~restriction_satisfied
                        continue    
                    end
                end
                % if all restriction passed, keep the draw and record
                H = self.mcmc_chol_Sigma(:,:,j) * Q;
                inv_H = Q' * la.invert_lower_triangular_matrix(self.mcmc_chol_Sigma(:,:,j));
                mcmc_H(:,:,i) = H;
                mcmc_inv_H(:,:,i) = inv_H;
                svar_index(i) = j;                
                if self.verbose
                    cu.progress_bar(i, self.iterations, 'Structural identification:');
                end
                i = i + 1;
            self.mcmc_H = mcmc_H;
            self.mcmc_inv_H = mcmc_inv_H;
            self.mcmc_Gamma = ones(self.iterations,self.n);
            self.svar_index = svar_index; 
            end
        end
        

        function svar_estimates(self)
            
            H_estimates = quantile(self.mcmc_H,0.5,3);
            Gamma_estimates = quantile(self.mcmc_Gamma,0.5,1);
            self.H_estimates = H_estimates;
            self.Gamma_estimates = Gamma_estimates;        
        end        
        
        
        function steady_state(self)
            
            ss = zeros(self.T,self.n,self.iterations);
            for j=1:self.iterations
                ss(:,:,j) = vu.steady_state(self.Z, self.mcmc_beta(:,:,j), self.n, self.m, self.p, self.T);
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Steady-state:');
                end
            end
            ss_estimates = vu.posterior_estimates(ss, self.credibility_level);
            self.steady_state_estimates = ss_estimates;
        end
        
       
        function fitted_and_residual(self)
        
            fitted = zeros(self.T,self.n,self.iterations);
            residual = zeros(self.T,self.n,self.iterations);
            for j=1:self.iterations
                [residual(:,:,j) fitted(:,:,j)] = vu.fit_and_residuals(self.Y, self.X, self.mcmc_beta(:,:,j));
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Fitted and residual:');
                end
            end
            fitted_estimates = vu.posterior_estimates(fitted, self.credibility_level);
            residual_estimates = vu.posterior_estimates(residual, self.credibility_level); 
            self.fitted_estimates = fitted_estimates;
            self.residual_estimates = residual_estimates;
            if self.structural_identification ~= 1
                structural_shocks = zeros(self.T,self.n,self.iterations);
                for j=1:self.iterations
                    index = self.svar_index(j);
                    structural_shocks(:,:,j) = vu.structural_shocks(residual(:,:,index), self.mcmc_inv_H(:,:,j));
                    if self.verbose
                        cu.progress_bar(j, self.iterations, 'Structural shocks:');
                    end
                end
                structural_shock_estimates = vu.posterior_estimates(structural_shocks, self.credibility_level);
                self.mcmc_structural_shocks = structural_shocks;
                self.structural_shock_estimates = structural_shock_estimates;
            end
        end
        
        
        function insample_criteria(self)
        
            insample_evaluation = vu.insample_evaluation_criteria(self.Y, ...
                                  self.residual_estimates(:,:,1), self.T, self.k);
            if self.verbose
                cu.progress_bar_complete('In-sample evaluation criteria:');
            end
            self.insample_evaluation = insample_evaluation;
        end
        
        
        function make_forecast(self, h, Z_p)         
        
            % make regressors
            [Z_p Y] = vu.make_forecast_regressors(Z_p, self.Y, h, self.p, self.T,...
                       self.exogenous, self.constant, self.trend, self.quadratic_trend);
            % initiate storage and loop over iterations
            mcmc_forecast = zeros(h,self.n,self.iterations);
            for j=1:self.iterations
                % make MCMC simulation for beta and Sigma
                mcmc_forecast(:,:,j) = vu.forecast(self.mcmc_beta(:,:,j),...
                             self.mcmc_chol_Sigma(:,:,j), h, Z_p, Y, self.n);
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Forecasts:');
                end
            end
            self.mcmc_forecast = mcmc_forecast;
        end
            
            
        function forecast_posterior_estimates(self, credibility_level)           
            
            % obtain posterior estimates
            mcmc_forecast = self.mcmc_forecast;
            forecast_estimates = vu.posterior_estimates(mcmc_forecast, credibility_level);
            self.forecast_estimates = forecast_estimates;
        end

        
        function make_impulse_response_function(self, h)        

            mcmc_irf = zeros(self.n, self.n, h, self.iterations);
            for i=1:self.iterations
                % get regular impulse response function
                mcmc_irf(:,:,:,i) = vu.impulse_response_function(self.mcmc_beta(:,:,i), self.n, self.p, h);
                if self.verbose   
                    cu.progress_bar(i, self.iterations, 'Impulse response function:');
                end
            end
            self.mcmc_irf = mcmc_irf;
        end
        

        function make_exogenous_impulse_response_function(self, h)

            % create exogenous IRFs only if there are exogenous (other than constant, trend, quadratic trend)
            if ~isempty(self.exogenous)
                % identify the number of exogenous (other than constant, trend, quadratic trend)
                r = size(self.exogenous,2);
                mcmc_irf_exo = zeros(self.n, r, h, self.iterations);
                for i=1:self.iterations
                    % get exogenous IRFs
                    mcmc_irf_exo(:,:,:,i) = vu.exogenous_impulse_response_function(self.mcmc_beta(:,:,i), self.n, self.m, r, self.p, h);
                    if self.verbose
                        cu.progress_bar(i, self.iterations, 'Exogenous impulse response function:');
                    end
                end
            else
                mcmc_irf_exo = [];
            end
            self.mcmc_irf_exo = mcmc_irf_exo;
        end

        
        function make_structural_impulse_response_function(self, h)

            if self.structural_identification == 1
                if self.verbose
                    cu.progress_bar_complete('Structural impulse response function:');
                end
                self.mcmc_structural_irf = [];                
            else
                mcmc_structural_irf = self.mcmc_irf;
                for i=1:self.iterations
                    % get index in case svar comes from restrictions
                    index = self.svar_index(i);
                    % recover structural impulse response function
                    mcmc_structural_irf(:,:,:,i) = vu.structural_impulse_response_function(...
                                                   self.mcmc_irf(:,:,:,index),self.mcmc_H(:,:,i), self.n);
                    if self.verbose   
                        cu.progress_bar(i, self.iterations, 'Structural impulse response function:');
                    end
                end
                self.mcmc_structural_irf = mcmc_structural_irf;
            end
        end
        
        
        function irf_posterior_estimates(self, credibility_level)

            if self.structural_identification == 1
                mcmc_irf = self.mcmc_irf;
            else
                mcmc_irf = self.mcmc_structural_irf;
            end
            irf_estimates = vu.posterior_estimates_3d(mcmc_irf, credibility_level);
            if ~isempty(self.exogenous)
                exo_irf_estimates = vu.posterior_estimates_3d(self.mcmc_irf_exo, credibility_level);
            else
                exo_irf_estimates = [];
            end
            self.irf_estimates = irf_estimates;
            self.exo_irf_estimates = exo_irf_estimates;
        end
        

        function make_forecast_error_variance_decomposition(self, h)

            % if no structural identification, FEVD is not computed
            if self.structural_identification == 1
                if self.verbose
                    cu.progress_bar_complete('Forecast error variance decomposition:');
                end
                self.mcmc_fevd = [];                
            % if there is some structural identification, proceed
            else
                % recover H from MCMC record, and Gamma if not identity
                if self.structural_identification == 3
                    mcmc_Gamma = num2cell(self.mcmc_Gamma,2);
                else
                    mcmc_Gamma = cell(self.iterations,1);
                end
                % initiate MCMC storage for FEVD and loop over iterations
                mcmc_fevd = zeros(self.n, self.n, h, self.iterations);                
                has_irf = ~isempty(self.mcmc_structural_irf) && size(self.mcmc_structural_irf, 3) >= h;
                for i=1:self.iterations
                    % recover structural IRF or estimate them
                    if has_irf
                        structural_irf = self.mcmc_structural_irf(:,:,1:h,i);
                    else
                        index = self.svar_index(i);
                        irf = vu.impulse_response_function(self.mcmc_beta(:,:,index), self.n, self.p, h);
                        structural_irf = vu.structural_impulse_response_function(irf, self.mcmc_H(:,:,i), self.n);
                    end
                    % recover fevd
                    mcmc_fevd(:,:,:,i) = vu.forecast_error_variance_decomposition(structural_irf, mcmc_Gamma{i}, self.n, h); 
                    if self.verbose   
                        cu.progress_bar(i, self.iterations, 'Forecast error variance decomposition:');
                    end
                end
                self.mcmc_fevd = mcmc_fevd;
            end
        end
        
        
        function fevd_posterior_estimates(self, credibility_level)

            if self.structural_identification == 1
                self.fevd_estimates = [];
            else
                mcmc_fevd = self.mcmc_fevd;
                fevd_estimates = vu.posterior_estimates_3d(mcmc_fevd, credibility_level);
                normalized_fevd_estimates = vu.normalize_fevd_estimates(fevd_estimates);
                self.fevd_estimates = normalized_fevd_estimates;
            end
        end

        
        function make_historical_decomposition(self)

            % if no structural identification, historical decomposition is not computed
            if self.structural_identification == 1
                if self.verbose
                    cu.progress_bar_complete('Historical decomposition:');
                end
                self.mcmc_hd = [];     
            % if there is some structural identification, proceed
            else
                % initiate MCMC storage for HD and loop over iterations
                mcmc_hd = zeros(self.n, self.n, self.T, self.iterations); 
                has_irf = ~isempty(self.mcmc_structural_irf) && size(self.mcmc_structural_irf, 3) >= self.T;
                has_structural_shocks = ~isempty(self.mcmc_structural_shocks);
                % loop over iterations
                for i=1:self.iterations
                    index = self.svar_index(i);
                    % recover structural IRF or estimate them
                    if has_irf
                        structural_irf = self.mcmc_structural_irf(:,:,1:self.T,i);
                    else
                        irf = vu.impulse_response_function(self.mcmc_beta(:,:,index), self.n, self.p, self.T);
                        structural_irf = vu.structural_impulse_response_function(irf, self.mcmc_H(:,:,i), self.n);
                    end
                    % recover structural shocs or estimate them
                    if has_structural_shocks
                        structural_shocks = self.mcmc_structural_shocks(:,:,i); 
                    else
                        [E ~] = vu.fit_and_residuals(self.Y, self.X, self.mcmc_beta(:,:,index));
                        structural_shocks = vu.structural_shocks(E, self.mcmc_inv_H(:,:,i));
                    end
                    % get historical decomposition
                    mcmc_hd(:,:,:,i) = vu.historical_decomposition(structural_irf, structural_shocks, self.n, self.T);          
                    if self.verbose
                        cu.progress_bar(i, self.iterations, 'Historical decomposition:');
                    end
                end
                self.mcmc_hd = mcmc_hd;
            end
        end        
        
        
        function hd_posterior_estimates(self, credibility_level)

            if self.structural_identification == 1
                self.hd_estimates = [];
            else
                mcmc_hd = self.mcmc_hd;
                hd_estimates = vu.posterior_estimates_3d(mcmc_hd, credibility_level);
                self.hd_estimates = hd_estimates;
            end
        end      
        
        
        function make_conditional_forecast(self, h, conditions, Z_p)        

            % make regressors Z_p and Y
            [Z_p Y] = vu.make_forecast_regressors(Z_p, self.Y, h, self.p, self.T,...
                      self.exogenous, self.constant, self.trend, self.quadratic_trend);
            % make conditional forecast regressors y_bar, Z, omega and gamma_00
            [y_bar Q omega gamma_00] = vu.conditional_forecast_regressors_1(conditions, h, Y, self.n, self.p);
            % initiate storage and loop over iterations
            mcmc_conditional_forecast = zeros(h,self.n,self.iterations);
            for i=1:self.iterations
                % recover iteration-specific regressors
                [mu F K Upsilon_00] = vu.conditional_forecast_regressors_2(self.mcmc_beta(:,:,i), ...
                                      self.mcmc_Sigma(:,:,i), conditions, Z_p, self.n, self.m, self.p, h);
                % run Carter Kohn algorithm to obtain conditional forecasts
                bss = BayesianStateSpaceSampler(y_bar, Q, omega, mu, F, K, 'z_00', gamma_00, ...
                      'Upsilon_00', Upsilon_00, 'kalman_type', 'conditional_forecast');
                bss.carter_kohn_algorithm();
                mcmc_conditional_forecast(:,:,i) = bss.Z(:,1:self.n);
                if self.verbose
                    cu.progress_bar(i, self.iterations, 'Conditional forecasts:');
                end
            self.mcmc_conditional_forecast = mcmc_conditional_forecast;
            end
        end
                
        
        function [shock_type] = check_shock_type(self, h, conditions, shocks)

            % check for structural identification
            if self.structural_identification == 1 
                if self.verbose
                    cu.progress_bar_complete('Conditional forecasts:');
                end
                shock_type = 'none';
            else
                % identify shocks
                if sum(shocks) == self.n
                    shock_type = 'all_shocks';
                else
                    shock_type = 'shock-specific';  
                end
            end
        end
        
        
        function make_structural_conditional_forecast(self, h, conditions, shocks, Z_p, shock_type)        

            % if there is an issue, return empty mcmc matrix
            if isequal(shock_type, 'none')
                self.mcmc_conditional_forecast = [];
                if self.verbose
                    cu.progress_bar_complete('Conditional forecasts:');
                end
            % if condition type is well defined, proceed
            else
                % make regressors Z_p and Y
                [Z_p Y] = vu.make_forecast_regressors(Z_p, self.Y, h, self.p, self.T,...
                          self.exogenous, self.constant, self.trend, self.quadratic_trend);
                has_irf = ~isempty(self.mcmc_structural_irf) && size(self.mcmc_structural_irf, 3) >= h;
                % make conditional forecast regressors R, y_bar and omega
                [R y_bar omega] = vu.conditional_forecast_regressors_3(conditions, h, self.n);
                if isequal(shock_type, 'shock-specific')
                    [P non_generating] = vu.conditional_forecast_regressors_5(shocks, h, self.n);
                end
                % initiate storage and loop over iterations
                mcmc_conditional_forecast = zeros(h,self.n,self.iterations);
                for i=1:self.iterations 
                    index = self.svar_index(i);
                    % make predictions, absent shocks
                    f = vu.linear_forecast(self.mcmc_beta(:,:,index), h, Z_p, Y, self.n);
                    % recover structural IRF or estimate them
                    if has_irf
                        structural_irf = self.mcmc_structural_irf(:,:,1:h,i);
                    else
                        irf = vu.impulse_response_function(self.mcmc_beta(:,:,index), self.n, self.p, h);
                        structural_irf = vu.structural_impulse_response_function(irf, self.mcmc_H(:,:,i), self.n);
                    end
                    % recover iteration-specific regressors
                    M = vu.conditional_forecast_regressors_4(structural_irf, self.n, h);
                    % get posterior mean and variance, depending on condition type    
                    if isequal(shock_type, 'all_shocks')
                        [mu_hat Omega_hat] = vu.conditional_forecast_posterior(y_bar, f, M, R, self.mcmc_Gamma(i,:)', omega, self.n, h);
                    elseif isequal(shock_type, 'shock-specific')
                        Gamma_nd = vu.conditional_forecast_regressors_6(self.mcmc_Gamma(i,:)', non_generating, h);                      
                        [mu_hat Omega_hat] = vu.shock_specific_conditional_forecast_posterior(...
                                             y_bar, f, M, R, P, self.mcmc_Gamma(i,:)', Gamma_nd, omega, self.n, h);  
                    end
                    % sample values
                    mcmc_conditional_forecast(:,:,i) = reshape(rn.multivariate_normal(mu_hat, Omega_hat)',[self.n,h])';
                    if self.verbose
                        cu.progress_bar(i, self.iterations, 'Conditional forecasts:');
                    end
                end
                self.mcmc_conditional_forecast = mcmc_structural_conditional_forecast;
            end
        end  

        
        function conditional_forecast_posterior_estimates(self, credibility_level)

            if isempty(self.mcmc_conditional_forecast)
                self.conditional_forecast_estimates = [];
            else
                mcmc_conditional_forecast = self.mcmc_conditional_forecast;
                conditional_forecast_estimates = vu.posterior_estimates(mcmc_conditional_forecast, credibility_level);
                self.conditional_forecast_estimates = conditional_forecast_estimates; 
            end
        end           


    end
    
end
