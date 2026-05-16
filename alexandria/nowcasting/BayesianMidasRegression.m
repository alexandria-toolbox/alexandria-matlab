classdef BayesianMidasRegression < handle


    % Bayesian MIDAS regression, developed in chapter 19
    % 
    % Parameters:
    % -----------
    % endogenous : matrix of size (n_obs,1)
    %     endogenous variables, defined in (6.19.1)
    % 
    % exogenous : matrix of size (n_obs,n_exogenous)
    %     exogenous variables, defined in (6.19.1)       
    % 
    % endogenous_lags : int, default = 1
    %     number of endogenous lags, defined in (6.19.1)
    % 
    % exogenous_lags : int or matrix of size (n_exogenous,1), default = 4
    %     number of exogenous lags, defined in (6.19.1) 
    % 
    % representation : char, default = 'unrestricted'
    %     applicable model representation, among 'unrestricted', 'almon', 'fourier' 
    %
    % prior_type : char, default = 'minnesota'
    %     applicable prior, among 'minnesota', 'horseshoe', 'lasso', 'almon', 'fourier'     
    % 
    % omega1 : float, default = 0.01
    %     overall tightness hyperparameter on alpha, defined in (6.19.10)
    % 
    % omega2 : float, default = 1
    %     cross-variable shrinkage on alpha, defined in (6.19.10)
    % 
    % upsilon1 : float, default = 0.1
    %     overall tightness hyperparameter on beta and delta, defined in (6.19.11)
    % 
    % upsilon2 : float, default = 1
    %     cross-variable shrinkage on beta and delta, defined in (6.19.11)
    %
    % polynomial_order : int, default = 3
    %     order of Alman or fourier polynomial, defined in (6.19.65) and (6.19.71)     
    % 
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    % 
    % burnin : int, default = 1000
    %     number of Gibbs sampler burn-in replications  
    % 
    % iterations : int, default = 2000
    %     number of Gibbs sampler replications   
    % 
    % verbose : bool, default = False
    %     if True, displays a progress bar 
    %
    %
    % Attributes
    % ----------
    % endogenous : matrix of size (n_obs,1)
    %     endogenous variables, defined in (6.19.1)
    % 
    % exogenous : matrix of size (n_obs,n_exogenous)
    %     exogenous variables, defined in (6.19.1)       
    % 
    % endogenous_lags : int, default = 1
    %     number of endogenous lags, defined in (6.19.1)
    % 
    % exogenous_lags : int or matrix of size (n_exogenous,1), default = 4
    %     number of exogenous lags, defined in (6.19.1) 
    % 
    % representation : char, default = 'unrestricted'
    %     applicable model representation, among 'unrestricted', 'almon', 'fourier' 
    %
    % prior_type : char, default = 'minnesota'
    %     applicable prior, among 'minnesota', 'horseshoe', 'lasso', 'almon', 'fourier'     
    % 
    % omega1 : float, default = 0.01
    %     overall tightness hyperparameter on alpha, defined in (6.19.10)
    % 
    % omega2 : float, default = 1
    %     cross-variable shrinkage on alpha, defined in (6.19.10)
    % 
    % upsilon1 : float, default = 0.1
    %     overall tightness hyperparameter on beta and delta, defined in (6.19.11)
    % 
    % upsilon2 : float, default = 1
    %     cross-variable shrinkage on beta and delta, defined in (6.19.11)
    %
    % polynomial_order : int, default = 3
    %     order of Alman or fourier polynomial, defined in (6.19.65) and (6.19.71)     
    % 
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    % 
    % burnin : int, default = 1000
    %     number of Gibbs sampler burn-in replications  
    % 
    % iterations : int, default = 2000
    %     number of Gibbs sampler replications   
    % 
    % verbose : bool, default = False
    %     if True, displays a progress bar 
    % 
    % p : int
    %     number of endogenous lags, defined in (6.19.1)
    % 
    % y : matrix of size (T,1)
    %     endogenous variables, defined in (6.19.7)        
    % 
    % Y : matrix of size (T,p)
    %     lagged endogenous variables, defined in (6.19.7)          
    % 
    % T : int
    %     number of realigned sample periods, defined in (6.19.1)        
    % 
    % n : int
    %     number of exogenous variables, defined in (6.19.1)    
    % 
    % p_ : matrix of size (n,1)
    %     number of exogenous lags p_i, defined in (6.19.1)     
    % 
    % X_ : cell of size (n)
    %     cell of endogenous regressors X_i, defined in (6.19.7)
    % 
    % X : matrix of size (T,l)
    %     full matrix of regressors X, defined in (6.19.7)
    %     
    % l : int
    %     number of regression coefficients, defined in (6.19.5)  
    % 
    % Xp : matrix of size (l,1)
    %     matrix of exogenous predictors, defined in (6.19.75)         
    % 
    % q_ : matrix of size (n,1)
    %     polynomial orders q_i, defined in (6.19.65) and (6.19.71)        
    % 
    % Q_ : cell of size (n)
    %     cell of parsimonious representation matrices Q_i, defined in (6.19.76) and (6.19.72)         
    % 
    % Z_ : cell of size (n)
    %     cell of parsimonious regressors Z_i, defined in (6.19.70)  
    % 
    % Z : matrix of size (T,n_parsimonious)
    %     full matrix of parsimonious regressors Z, defined in (6.19.70)     
    % 
    % kappa : float
    %     prior shape on sigma, defined in (6.19.14)
    % 
    % lamda : float
    %     prior scale on sigma, defined in (6.19.14)
    % 
    % f : float
    %     prior shape on mu, defined in (6.19.50)
    % 
    % g : float
    %     prior scale on mu, defined in (6.19.50)
    % 
    % mcmc_beta : cell of size (f_periods)
    %     cell of mcmc values for regression coefficients beta_i, for each forecast period
    % 
    % mcmc_sigma : cell of size (f_periods)
    %     cell of mcmc values for sigma, for each forecast period
    % 
    % mcmc_forecast : cell of size (f_periods)
    %     cell of mcmc values for predictions, for each forecast period
    % 
    % forecast_estimates : cell of size (f_periods)
    %     cell of forecast estimates, for each forecast period
    %     entries are median, lower bounds, upper bounds
    % 
    % beta_estimates: cell of size (n)
    %     median, lower bound, upper bound and variance of regression coefficeints beta_i
    % 
    % sigma_estimates: matrix of size (4,1)
    %     median, lower bound, upper bound and variance of residual variance sigma
    % 
    % fitted_estimates: matrix of size (T,3)
    %     estimates of in-sample fit
    % 
    % residual_estimates: matrix of size (T,3)
    %     estimates of in-sample residuals
    % 
    % insample_evaluation: struct of size(3)
    %     estimates of in-sample evaluation criteria
    % 
    % forecast_evaluation_criteria: struct
    %     structure of forecast evaluation criteria  
    % 
    % 
    % Methods
    % ----------
    % estimate
    % insample_fit
    % forecast 
    % forecast_evaluation


    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
     
    
    properties (GetAccess = public, SetAccess= protected)
        endogenous
        exogenous
        endogenous_lags
        exogenous_lags
        representation
        prior_type
        omega1
        omega2
        upsilon1
        upsilon2        
        polynomial_order
        credibility_level
        iterations
        burnin
        verbose
        p
        y
        Y
        T
        n
        p_
        X_
        X
        l
        Xp
        q_
        Q_
        Z_
        Z
        kappa
        lamda
        f
        g
        mcmc_beta
        mcmc_sigma
        mcmc_forecast      
        forecast_estimates
        beta_estimates
        sigma_estimates
        fitted_estimates
        residual_estimates
        insample_evaluation
        forecast_evaluation_criteria
    end
    
    
    properties (GetAccess = private, SetAccess = private)
        non_nan_entries
        one_T
        ll
        inv_V
        endogenous_indices
    end
    

    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = BayesianMidasRegression(endogenous, exogenous, varargin)
            
            % constructor for the BayesianMidasRegression class
            
            % allow for optional arguments
            parser = inputParser;
            default_endogenous_lags = 1;
            default_exogenous_lags = 4;
            default_representation = 'unrestricted';
            default_prior_type = 'minnesota';
            default_omega1 = 0.01;
            default_omega2 = 1;
            default_upsilon1 = 0.1;
            default_upsilon2 = 1;            
            default_polynomial_order = 2;
            default_credibility_level = 0.95;
            default_iterations = 2000;
            default_burnin = 1000;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addRequired(parser, 'exogenous');
            addParameter(parser, 'endogenous_lags', default_endogenous_lags);
            addParameter(parser, 'exogenous_lags', default_exogenous_lags);
            addParameter(parser, 'representation', default_representation);
            addParameter(parser, 'prior_type', default_prior_type);
            addParameter(parser, 'omega1', default_omega1); 
            addParameter(parser, 'omega2', default_omega2); 
            addParameter(parser, 'upsilon1', default_upsilon1); 
            addParameter(parser, 'upsilon2', default_upsilon2);             
            addParameter(parser, 'polynomial_order', default_polynomial_order);
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'burnin', default_burnin);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'verbose', default_verbose);
            parse(parser, endogenous, exogenous, varargin{:});
            self.endogenous = endogenous;
            self.exogenous = exogenous;
            self.endogenous_lags = parser.Results.endogenous_lags;
            self.exogenous_lags = parser.Results.exogenous_lags;
            self.representation = parser.Results.representation;
            self.prior_type = parser.Results.prior_type;
            self.omega1 = parser.Results.omega1;
            self.omega2 = parser.Results.omega2;
            self.upsilon1 = parser.Results.upsilon1;
            self.upsilon2 = parser.Results.upsilon2;            
            self.polynomial_order = parser.Results.polynomial_order;
            self.credibility_level = parser.Results.credibility_level;
            self.burnin = parser.Results.burnin;
            self.iterations = parser.Results.iterations;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
            % define prior values
            self.prior();     
            % initialize record elements
            self.initialize_records();         
        end


        function estimate(self)

            % estimate()
            % generates model estimates and nowcasts for the Bayesian Midas model
            % 
            % parameters:
            % none
            % 
            % returns:
            % none
            
            % update recording parameters
            self.update_records(1);
            % run MCMC algorithm for regression parameters
            self.parameter_mcmc(1);
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
        end


        function insample_fit(self)
        
            % insample_fit()
            % generates in-sample fit and residuals
            % 
            % parameters:
            % none
            % 
            % returns:
            % none          
            
            % compute fitted and residuals
            self.fitted_and_residual();
            % compute in-sample criteria
            self.insample_criteria();
        end


        function forecast(self, h, credibility_level)
            
            % forecast(h, credibility_level)
            % generates model estimates and nowcasts for the Bayesian Midas model
            % 
            % parameters:
            % h : int
            %     number of forecast periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % none       
    
            % update recording parameters
            self.update_records(h);
            % run MCMC algorithm for regression parameters if not already done
            if size(self.mcmc_beta{h},1) == 0
                self.parameter_mcmc(h);
            end
            % run mcmc algorithm for predictions if not already done
            if size(self.mcmc_forecast{h},1) == 0
                self.make_forecast(h);
                self.forecast_posterior_estimates(h, credibility_level);
            end
        end


        function forecast_evaluation(self, y)
            
            % forecast_evaluation(y)
            % forecast evaluation criteria for the linear regression model
            % 
            % parameters:
            % y : ndarray of shape (m,)
            %     array of realised values for forecast evaluation
            % 
            % returns:
            % forecast_evaluation_criteria: dict
            %     dictionary of forecast evaluation criteria
            
            % recover forecasts
            iterations = self.iterations;
            indices = find(~cellfun(@isempty,self.mcmc_forecast));
            m = size(indices,1);
            y_hat = zeros(m,1);
            y_ = zeros(m,1);
            mcmc_forecast = zeros(m,iterations);
            for i=1:m
                y_hat(i) = self.forecast_estimates{indices(i)}(1);
                y_(i) = y(indices(i));
                mcmc_forecast(i,:) = self.mcmc_forecast{indices(i)};
            end
            % get regular forecast evaluation criteria
            standard_evaluation_criteria = ru.forecast_evaluation_criteria(y_hat, y_); 
            % obtain Bayesian forecast evaluation criteria
            bayesian_evaluation_criteria = self.bayesian_forecast_evaluation_criteria(y, mcmc_forecast, iterations, m);  
            % merge structures
            forecast_evaluation_criteria = iu.concatenate_structures(standard_evaluation_criteria, bayesian_evaluation_criteria);                   
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;
        end


    end


    methods (Access = protected, Hidden = true)


        function make_regressors(self)
                
            % generates regressors using frequency alignment

            % determine initial observation
            self.make_initial_observation();
            % make endogenous regressors
            self.make_endogenous_regressors();
            % make exogenous regressors
            self.make_exogenous_regressors();
            % make predictors
            self.make_predictors();
            % make parsimonious representation
            self.make_parsimonious_representation();
        end


        function make_initial_observation(self)

            % initial endogenous variable observation so that lags are well defined

            self.n = size(self.exogenous,2);
            if isscalar(self.exogenous_lags)
                self.p_ = self.exogenous_lags * ones(self.n,1);
            else
                self.p_ = self.exogenous_lags;
            end
            indices = zeros(self.n+1,1);
            for i=1:self.n
                [non_nan_endogenous non_nan_entries] = la.dropna(self.exogenous(:,i));
                indices(i) = non_nan_entries(self.p_(i)+1);  
            end
            self.p = self.endogenous_lags;
            [non_nan_endogenous non_nan_entries] = la.dropna(self.endogenous);
            indices(end) = non_nan_entries(self.p+1);
            self.endogenous_indices = non_nan_entries(non_nan_entries >= max(indices));
        end


        function make_endogenous_regressors(self)
            
            % endogenous regressors y and Y defined in (6.19.3)

            self.T = size(self.endogenous_indices,1);
            self.y = zeros(self.T,1);
            self.Y = zeros(self.T,self.p);
            for t=1:self.T
                [y_t ~] = la.dropna(self.endogenous(1:self.endogenous_indices(t)));
                self.y(t) = y_t(end);
                if self.p > 0
                    self.Y(t,:) = flipud(y_t(end-self.p:end-1));
                end
            end
        end


        function make_exogenous_regressors(self)
            
            % endogenous regressors 1_T and X_i defined in (6.19.3)
            
            self.one_T = ones(self.T,1);
            X = [self.one_T, self.Y];
            X_ = cell(self.n,1);
            for i=1:self.n
                temp = self.exogenous(:,i);
                X_i = zeros(self.T,self.p_(i)+1);
                for t=1:self.T
                    [X_it ~] = la.dropna(temp(1:self.endogenous_indices(t)));
                    X_it = flipud(X_it(end-self.p_(i):end));
                    X_i(t,:) = X_it;
                end
                X_{i} = X_i;
                X = [X X_i];
            end
            self.X_ = X_;
            self.X = X;
            self.l = 1 + self.p + sum(self.p_) + self.n;
        end


        function make_predictors(self)
            
            % endogenous and exogenous predictors

            Xp_ = cell(self.n,1);
            for i=1:self.n
                [temp ~] = la.dropna(self.exogenous(:,i));
                Xp_i = temp(end-self.p_(i):end);
                Xp_{i} = Xp_i;
            end
            if self.p == 0
                Yp = zeros(0,1);
            else
                [temp ~] = la.dropna(self.endogenous);
                Yp = flipud(temp(end-self.p+1:end));
            end
            self.Xp = [1;Yp;cell2mat(Xp_)];
        end
            

        function make_parsimonious_representation(self)
            
        % parsimonious representations, defined in (6.19.6)-(6.19.11)

            if isscalar(self.polynomial_order)
                self.q_ = self.polynomial_order * ones(self.n,1);
            else
                self.q_ = self.polynomial_order;
            end
            Q_ = cell(self.n,1);
            Z_ = cell(self.n,1);
            Z = [self.one_T self.Y];
            for i=1:self.n
                if isequal(self.representation,'almon')
                    temp_1 = (0:self.p_(i))' .* ones(1,self.q_(i)+1);
                    temp_2 = ones(self.p_(i)+1,1) * (0:self.q_(i));
                    Q_i = temp_1 .^ temp_2;
                elseif isequal(self.representation,'fourier')
                    temp_1 = ones(self.p_(i)+1,1);
                    temp_2 = 2 * pi / (self.p_(i) + 1) * (0:self.p_(i))' * (1:self.q_(i));
                    temp_3 = cos(temp_2);
                    temp_4 = sin(temp_2);
                    Q_i = [temp_1 temp_3 temp_4];
                end
                if isequal(self.representation,'almon') || isequal(self.representation,'fourier')
                    Q_{i} = Q_i;
                    Z_i = self.X_{i} * Q_i;
                    Z_{i} = Z_i;
                    Z = [Z Z_i];
                end
            end
            self.Q_ = Q_;
            self.Z_ = Z_;
            if isequal(self.representation,'unrestricted')
                self.Z = [];
                self.ll = 1 + self.p + self.n + sum(self.p_);
            elseif isequal(self.representation,'almon')
                self.Z = Z;
                self.ll = 1 + self.p + self.n + sum(self.q_);
            elseif isequal(self.representation,'fourier')
                self.Z = Z;
                self.ll = 1 + self.p + self.n + 2 * sum(self.q_);
            end
        end
     
            
        function prior(self)
            
            % creates prior elements
            
            % Minnesota prior hyperparameters, defined in (6.19.12)
            if isequal(self.prior_type, 'minnesota')
                u = 10000;
                w = (self.omega1 ./ (1:self.p) .^ self.omega2) .^ 2;
                v = [u w];
                for i=1:self.n
                    if isequal(self.representation, 'unrestricted')
                        v_i = (self.upsilon1 ./ (0:self.p_(i)) .^ self.upsilon2) .^ 2;
                    elseif isequal(self.representation, 'almon')
                        v_i = (self.upsilon1 ./ (0:self.q_(i)) .^ self.upsilon2) .^ 2;
                    elseif isequal(self.representation, 'fourier')
                        v_i = (self.upsilon1 ./ [(0:self.q_(i)) (1:self.q_(i))] .^ self.upsilon2) .^ 2;
                    end
                    v_i(1) = (2 * self.upsilon1) ^ 2;
                    v = [v v_i];
                end
                inv_V = diag(1 ./ v);
                self.inv_V = inv_V;
            % Bayesian lasso hyperparameters, defined in (6.19.50)    
            elseif isequal(self.prior_type, 'lasso')
                self.f = 0.001;
                self.g = 0.001;
            end
            % generate kappa and lambda, defined in (6.19.14), (6.19.27) and (6.19.51)
            self.kappa = 0.001;
            self.lamda = 0.001;
        end


        function initialize_records(self)
            
            % initialize recording elements

            self.mcmc_beta = {};
            self.mcmc_sigma = {};
            self.mcmc_forecast = {};      
            self.forecast_estimates = {};
        end


        function update_records(self, horizon)
                
            % update recording elements

            length = numel(self.mcmc_beta);
            if length < horizon
                self.mcmc_beta = resize(self.mcmc_beta, [horizon 1]);
                self.mcmc_sigma = resize(self.mcmc_sigma, [horizon 1]);
                self.mcmc_forecast = resize(self.mcmc_forecast, [horizon 1]);
                self.forecast_estimates = resize(self.forecast_estimates, [horizon 1]);
            end
        end
    

        function parameter_mcmc(self, horizon)
            
            % posterior distributions for parameters from algorithms 19.1-19.4
            
            % unpack
            representation = self.representation;
            prior_type = self.prior_type;
            n = self.n;
            l = self.l;
            ll = self.ll;
            p = self.p;
            p_ = self.p_;
            q_ = self.q_;
            Q_ = self.Q_;
            kappa = self.kappa;
            lamda = self.lamda;
            if isequal(prior_type, 'minnesota')
                inv_V = self.inv_V;
            elseif isequal(prior_type, 'lasso')
                f = self.f;
                g = self.g;
            end
            burnin = self.burnin;
            iterations = self.iterations;
            verbose = self.verbose;
            if verbose && horizon == 1
                verbose_string = 'Model parameters, 1 period ahead:';
            elseif verbose && horizon > 1
                verbose_string = ['Model parameters, ' num2str(horizon) ' periods ahead:'];
            end
            
            % initialize storage
            mcmc_delta = zeros(ll,iterations); 
            mcmc_sigma = zeros(iterations,1);

            % realign regressors
            [y X T] = self.realign_regressors(horizon);

            % prior, posterior and initial values, depending on prior type
            if isequal(prior_type, 'minnesota')
                [XX Xy kappa_bar sigma] = self.minnesota_parameters(y, X, T, kappa);
                
            elseif isequal(prior_type, 'horseshoe')
                [XX Xy gamma_nu gamma_tau gamma_eta gamma_psi kappa_bar sigma tau2 psi2] ...
                    = self.horseshoe_parameters(y, X, T, ll, kappa);
                
            elseif isequal(prior_type, 'lasso')
                [inv_XX b_bar f_bar kappa_bar sigma xi mu] = ...
                    self.lasso_parameters(y, X, T, ll, kappa, f);
            end

            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)  

                % move on with Minnesota prior if selected
                if isequal(prior_type, 'minnesota')

                    % step 2: draw beta_i
                    [beta] = self.draw_minnesota_beta(Xy, XX, sigma, inv_V);

                    % step 3: draw sigma
                    [sigma] = self.draw_minnesota_sigma(kappa_bar, y, X, beta, lamda);

                % else, move on with horseshoe prior if selected
                elseif isequal(prior_type, 'horseshoe')
    
                    % step 4: draw beta
                    [beta beta2] = self.draw_horseshoe_beta(Xy, XX, tau2, psi2, sigma);

                    % step 5: draw nu
                    [nu] = self.draw_nu(gamma_nu, tau2);

                    % step 6: draw tau2
                    [tau2] = self.draw_tau2(gamma_tau, nu, sigma, beta2, psi2);
    
                    % step 7: draw eta
                    [eta] = self.draw_eta(gamma_eta, psi2, ll);

                    % step 8: draw psi2
                    [psi2] = self.draw_psi2(gamma_psi, eta, sigma, beta2, tau2, ll);

                    % step 9: draw sigma
                    [sigma] = self.draw_horseshoe_sigma(kappa_bar, y, X, beta, beta2, tau2, psi2, lamda);

                % else, move on with Bayesian lasso if selected
                elseif isequal(prior_type, 'lasso')

                    % step 4: draw beta
                    [beta beta2] = self.draw_lasso_beta(inv_XX, b_bar, sigma, xi);

                    % step 5: draw xi
                    [xi] = self.draw_xi(beta, mu, sigma, ll);

                    % step 6: draw mu
                    [mu] = self.draw_mu(f_bar, g, xi);
    
                    % step 7: draw sigma
                    [sigma] = self.draw_lasso_sigma(kappa_bar, y, X, beta, beta2, xi, lamda);
                end

                % save if burn is exceeded
                if iteration > burnin

                    % save parameter values
                    mcmc_delta(:,iteration-burnin) = beta;
                    mcmc_sigma(iteration-burnin) = sigma;                    
                end

                % if verbose, display progress bar
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, verbose_string); 
                end

                % update iterations
                iteration = iteration + 1;
            end

            % recover beta if parsimonious representation was selected
            mcmc_beta = self.recover_beta_from_parsimonious(mcmc_delta, ...
                        iterations, representation, n, l, p, p_, q_, Q_);

            % save as attributes
            self.mcmc_beta{horizon} = mcmc_beta;
            self.mcmc_sigma{horizon} = mcmc_sigma;
        end


        function [y X T] = realign_regressors(self, h)
            
            % realign regressors to forecast horizon
    
            y = self.y(h:end);
            T = size(y,1);
            one_T = self.one_T(1:T,:); 
            if self.p == 0
                Y = zeros(T,0);
            else
                Y = self.Y(1:T,:);
            end
            X = [one_T Y];
            for i=1:self.n
                if isequal(self.representation, 'unrestricted')
                    X_i = self.X_{i}(1:T,:);
                    X = [X X_i];
                elseif isequal(self.representation, 'almon') || isequal(self.representation, 'fourier')
                    Z_i = self.Z_{i}(1:T,:);
                    X = [X Z_i];
                end
            end
        end


        function [XX Xy kappa_bar sigma] = minnesota_parameters(self, y, X, T, kappa)
    
            % prior and posterior values for Minnesota
    
            XX = X' * X;
            Xy = X' * y;
            kappa_bar = (T + kappa) / 2;
            sigma = 1;
        end


        function [XX, Xy, gamma_nu, gamma_tau, gamma_eta, gamma_psi, kappa_bar, sigma, tau2, psi2_] ...
                = horseshoe_parameters(self, y, X, T, ll, kappa)
            
            % prior and posterior values for horeseshoe
            
            XX = X' * X;
            Xy = X' * y;
            gamma_nu = 1;
            gamma_tau = (ll + 1) / 2;
            gamma_eta = 1;
            gamma_psi = 1;
            kappa_bar = (T + ll + kappa) / 2;
            sigma = 1;
            tau2 = 5;
            psi2_ = 5 * ones(ll,1);    
        end


        function [inv_XX b_bar f_bar kappa_bar sigma xi mu] ...
                = lasso_parameters(self, y, X, T, ll, kappa, f)
            
            % prior and posterior values for lasso
            
            [inv_XX] = la.robust_covariance_matrix(X);
            Xy = X' * y;
            b_bar = inv_XX * Xy;
            f_bar = f + 2 * ll;
            kappa_bar = kappa + T + ll;
            sigma = 1;    
            xi = 5 * ones(ll,1);
            mu = 1;
        end


        function [beta] = draw_minnesota_beta(self, Xy, XX, sigma, inv_V)
            
            % draw beta from its conditional posterior defined in (6.19.16) 
            
            inv_V_bar = inv_V + XX / sigma;
            b_bar_temp = Xy / sigma;
            beta = rn.efficient_multivariate_normal(b_bar_temp, inv_V_bar);
        end


        function [sigma] = draw_minnesota_sigma(self, kappa_bar, y, X, beta, lamda)
            
            % draw sigma from its conditional posterior defined in (6.19.19)

            residuals = y - X * beta;
            lamda_bar = residuals' * residuals + lamda;
            sigma = rn.inverse_gamma(kappa_bar, lamda_bar);
        end


        function [beta beta2] = draw_horseshoe_beta(self, Xy, XX, tau2, psi2, sigma)
            
            % draw beta from its conditional posterior defined in (6.19.30)  

            inv_V_star = diag(1 ./ (psi2 * tau2)) + XX;
            V_star = la.invert_spd_matrix(inv_V_star);
            V_bar = sigma * V_star;
            b_bar = V_star * Xy;
            beta = rn.multivariate_normal(b_bar, V_bar);
            beta2 = beta .^ 2;
        end


        function [nu] = draw_nu(self, gamma_nu, tau2)
            
            % draw nu from its conditional posterior defined in (6.19.36) 

            phi_nu = 1 / tau2 + 1;
            nu = rn.inverse_gamma(gamma_nu, phi_nu);
        end


        function [tau2] = draw_tau2(self, gamma_tau, nu, sigma, beta2, psi2) 
     
            % draw tau2 from its conditional posterior defined in (6.19.33)
    
            phi_tau = 1 / nu + sum(beta2 ./ psi2) / (2 * sigma);
            tau2 = rn.inverse_gamma(gamma_tau, phi_tau) + 1e-12;
        end


        function [eta] = draw_eta(self, gamma_eta, psi2, ll)
            
            % draw eta from its conditional posterior defined in (6.19.42)

            eta = zeros(ll,1);
            phi_eta = 1 ./ psi2 + 1;
            for h=1:ll
                eta(h) = rn.inverse_gamma(gamma_eta, phi_eta(h));
            end
        end


        function [psi2] = draw_psi2(self, gamma_psi, eta, sigma, beta2, tau2, ll)
    
            % draw psi2 from its conditional posterior defined in (6.19.39)  

            psi2 = zeros(ll,1);
            phi_psi = 1 ./ eta + beta2 / (2 * sigma * tau2);
            for h=1:ll
                psi2(h) = rn.inverse_gamma(gamma_psi, phi_psi(h));
            end
        end


        function [sigma] = draw_horseshoe_sigma(self, kappa_bar, y, X, beta, beta2, tau2, psi2, lamda)
            
            % draw sigma from its conditional posterior defined in (6.19.45)

            residuals = y - X * beta;
            ratios = sum(beta2 ./ (tau2 * psi2));
            lamda_bar = (residuals' * residuals + ratios + lamda) / 2;
            sigma = rn.inverse_gamma(kappa_bar, lamda_bar);
        end


        function [beta beta2] = draw_lasso_beta(self, inv_XX, b_bar, sigma, xi)
            
            % draw beta from its conditional posterior defined in (6.19.54)       
      
            V_bar = sigma * inv_XX;
            chi = sqrt(sigma) * xi;
            beta = rn.truncated_multivariate_normal(b_bar, V_bar, -chi, chi);
            beta2 = beta.^2;
        end


        function [xi] = draw_xi(self, beta, mu, sigma, ll)
            
            % draw xi from its conditional posterior defined in (6.19.57)    

            sqrt_sigma = sqrt(sigma);
            xi = zeros(ll,1);
            for h=1:ll
                rho_h = abs(beta(h)) / sqrt_sigma;
                xi(h) = rn.gamma(1,mu) + rho_h;
            end
        end


        function [mu] = draw_mu(self, f_bar, g, xi)  
            
            % draw mu from its conditional posterior defined in (6.19.60)

            g_bar = g + sum(xi);
            mu = rn.inverse_gamma(f_bar,g_bar);
        end


        function [sigma] = draw_lasso_sigma(self, kappa_bar, y, X, beta, beta2, xi, lamda)
            
            % draw sigma from its conditional posterior defined in (6.19.63)

            residuals = y - X * beta;
            lamda_bar = (lamda + residuals' * residuals) / 2;
            iota_bar = max(beta2 ./ xi.^2);
            gamma_scale = 1 / lamda_bar;
            p_max = gammainc(1 / (iota_bar*gamma_scale), kappa_bar);
            sigma = 1 / (gamma_scale * gammaincinv(unifrnd(0,p_max), kappa_bar));
        end


        function [mcmc_beta] = recover_beta_from_parsimonious(self, mcmc_delta, ...
                               iterations, representation, n, l, p, p_, q_, Q_)
            
            % recover beta from delta using the equivalence defined in (6.19.66)
            
            if isequal(representation,'unrestricted')
                mcmc_beta = mcmc_delta;
            else
                mcmc_beta = zeros(l,iterations);
                mcmc_beta(1:1+p,:) = mcmc_delta(1:1+p,:);
                beta_index = cumsum([1+p;p_+1]);
                if isequal(representation,'almon')
                    delta_index = cumsum([1+p;q_+1]);
                elseif isequal(representation,'fourier')
                    delta_index = cumsum([1+p;2*q_+1]);
                end
                for i=1:n
                    Q_i = Q_{i};
                    mcmc_delta_i = mcmc_delta(delta_index(i)+1:delta_index(i+1),:);
                    mcmc_beta_i = Q_i * mcmc_delta_i;
                    mcmc_beta(beta_index(i)+1:beta_index(i+1),:) = mcmc_beta_i;
                end
            end
        end


        function parameter_estimates(self)
    
            % posterior estimates for midas regression parameters
            
            beta_estimates = zeros(self.l,4);
            beta_estimates(:,1) = quantile(self.mcmc_beta{1},0.5,2);
            beta_estimates(:,2) = quantile(self.mcmc_beta{1},(1-self.credibility_level)/2,2);
            beta_estimates(:,3) = quantile(self.mcmc_beta{1},(1+self.credibility_level)/2,2);    
            beta_estimates(:,4) = std(self.mcmc_beta{1},0,2);
            sigma_estimates = zeros(4,1);
            sigma_estimates(1) = quantile(self.mcmc_sigma{1},0.5);
            sigma_estimates(2) = quantile(self.mcmc_sigma{1},(1-self.credibility_level)/2);
            sigma_estimates(3) = quantile(self.mcmc_sigma{1},(1+self.credibility_level)/2);   
            sigma_estimates(4) = std(self.mcmc_sigma{1});
            self.beta_estimates = beta_estimates;
            self.sigma_estimates = sigma_estimates;
        end


        function fitted_and_residual(self)
            
            % fitted and residual values
            
            y = self.y;
            X = self.X;
            T = self.T;
            mcmc_beta = self.mcmc_beta{1};
            mcmc_fitted = zeros(T,self.iterations);
            mcmc_residual = zeros(T,self.iterations);
            for i=1:self.iterations
                y_hat = X * mcmc_beta(:,i);
                residuals = y - y_hat;
                mcmc_fitted(:,i) = y_hat;
                mcmc_residual(:,i) = residuals;
            end
            fitted_estimates = zeros(T,3);
            fitted_estimates(:,1) = quantile(mcmc_fitted,0.5,2);
            fitted_estimates(:,2) = quantile(mcmc_fitted,(1-self.credibility_level)/2,2);
            fitted_estimates(:,3) = quantile(mcmc_fitted,(1+self.credibility_level)/2,2);         
            residual_estimates = zeros(T,3);
            residual_estimates(:,1) = quantile(mcmc_residual,0.5,2);
            residual_estimates(:,2) = quantile(mcmc_residual,(1-self.credibility_level)/2,2);
            residual_estimates(:,3) = quantile(mcmc_residual,(1+self.credibility_level)/2,2);  
            self.fitted_estimates = fitted_estimates;
            self.residual_estimates = residual_estimates;
        end


        function insample_criteria(self)
            
            % in-sample fit evaluation criteria
    
            insample_evaluation = ru.insample_evaluation_criteria(self.y, self.residual_estimates(:,1), self.T, self.l);
            self.insample_evaluation = insample_evaluation;
        end


        function make_forecast(self, h) 
    
            % forecasts for the midas regression model       
    
            mcmc_forecast = zeros(self.iterations,1);
            if self.verbose && h == 1
                verbose_string = 'Forecasts, 1 period ahead:';
            elseif self.verbose && h > 1
                verbose_string = ['Forecasts, ' num2str(h) ' periods ahead:']; 
            end
            X = self.Xp;
            mcmc_beta = self.mcmc_beta{h};
            mcmc_sigma  = self.mcmc_sigma{h};
            for i=1:self.iterations
                beta = mcmc_beta(:,i);
                sigma  = mcmc_sigma(i);
                e = sqrt(sigma) * randn;
                yp = X' * beta + e;
                mcmc_forecast(i) = yp;
                if self.verbose
                    cu.progress_bar(i, self.iterations, verbose_string);
                end
            end
            self.mcmc_forecast{h} = mcmc_forecast;
        end


        function forecast_posterior_estimates(self, h, credibility_level)
        
            % posterior estimates for forecasts
    
            % obtain posterior estimates
            mcmc_forecast = self.mcmc_forecast{h};
            forecast_estimates = zeros(3,1);
            forecast_estimates(1) = quantile(mcmc_forecast, 0.5);
            forecast_estimates(2) = quantile(mcmc_forecast, (1-credibility_level)/2);
            forecast_estimates(3) = quantile(mcmc_forecast, (1+credibility_level)/2);
            self.forecast_estimates{h} = forecast_estimates;
        end


        function [bayesian_forecast_evaluation_criteria] = bayesian_forecast_evaluation_criteria(self, y, mcmc_forecast, iterations, m)
            
            % Bayesian forecast evaluation criteria from equations from equations (3.10.13) and (3.10.15)
    
            log_score = zeros(m,1);
            crps = zeros(m,1);
            for i = 1:m
                % get actual, prediction mean, prediction variance
                y_i = y(i);
                forecasts = mcmc_forecast(i,:);
                mu_i = mean(forecasts);
                sigma_i = var(forecasts);
                % get log score from equation (3.10.14)
                [log_pdf, ~] = su.normal_pdf(y_i, mu_i, sigma_i);
                log_score(i) = - log_pdf;
                % get CRPS from equation (3.10.17)
                term_1 = sum(abs(forecasts - y_i));
                term_2 = 0;
                for j = 1:iterations
                    term_2 = term_2 + sum(abs(forecasts(j) - forecasts));
                end
                crps(i) = term_1 / iterations - term_2 / (2 * iterations^2);
            end
            log_score = mean(log_score);
            crps = mean(crps);  
            bayesian_forecast_evaluation_criteria = struct;
            bayesian_forecast_evaluation_criteria.log_score = log_score;
            bayesian_forecast_evaluation_criteria.crps = crps;
        end


    end


end