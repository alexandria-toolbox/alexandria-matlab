classdef VectorAutoregressionResults < handle
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        regressors
        coefficient_index
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    


        function self = VectorAutoregressionResults()          
        end


    end
    
    
    methods (Access = protected, Hidden = true)
        

        function complete_var_information(self)
            % endogenous and exogenous variables
            if ~isfield(self.complementary_information, 'endogenous_variables')
                n_endo = size(self.model.endogenous,2);
                self.complementary_information.endogenous_variables = "y"+(1:n_endo);
            end
            if ~isfield(self.complementary_information, 'exogenous_variables')
                if isempty(self.model.exogenous)
                    self.complementary_information.exogenous_variables = "none";
                else
                    n_exo = size(self.model.exogenous,2);
                    self.complementary_information.exogenous_variables = "x"+(1:n_exo);
                end
            end
            % sample dates
            if ~isfield(self.complementary_information, 'dates')
                T = self.model.T;
                p = self.model.p;
                self.complementary_information.dates = (1-p:T)';
            end
            % forecast dates
            if ~isfield(self.complementary_information, 'forecast_dates')
                if ~isempty(self.model.forecast_estimates)
                    f_periods = size(self.model.forecast_estimates,1);
                    T = self.model.T;
                    self.complementary_information.forecast_dates = (T+1:T+f_periods)';
                else
                    self.complementary_information.forecast_dates = [];
                end
            end
            % conditional forecast dates
            if ~isfield(self.complementary_information, 'conditional_forecast_dates')
                if isfield(self.model, 'conditional_forecast_estimates') && ~isempty(self.model.conditional_forecast_estimates)
                    f_periods = size(self.model.conditional_forecast_estimates,1);
                    T = self.model.T;
                    self.complementary_information.conditional_forecast_dates = (T+1:T+f_periods)';
                else
                    self.complementary_information.conditional_forecast_dates = [];
                end
            end
            % proxy variables
            model_type = self.complementary_information.model_type;
            if model_type == 7 && ~isfield(self.complementary_information, 'proxy_variables')
                self.complementary_information.proxy_variables = "proxy_"+(1:size(self.model.proxys,2));
            end
            % VAR options
            if ~isfield(self.complementary_information, 'insample_fit')
                self.complementary_information.insample_fit = '—';
            end
            if ~isfield(self.complementary_information, 'constrained_coefficients')
                self.complementary_information.constrained_coefficients = '—';
            end
            if ~isfield(self.complementary_information, 'sums_of_coefficients')
                self.complementary_information.sums_of_coefficients = '—';
            end
            if ~isfield(self.complementary_information, 'initial_observation')
                self.complementary_information.initial_observation = '—';
            end
            if ~isfield(self.complementary_information, 'long_run')
                self.complementary_information.long_run = '—';
            end
            if ~isfield(self.complementary_information, 'stationary')
                self.complementary_information.stationary = '—';
            end
            if ~isfield(self.complementary_information, 'marginal_likelihood')
                self.complementary_information.marginal_likelihood = '—';
            end
            if ~isfield(self.complementary_information, 'hyperparameter_optimization')
                self.complementary_information.hyperparameter_optimization = '—';
            end
            if ~isfield(self.complementary_information, 'coefficients_file')
                self.complementary_information.coefficients_file = '—';
            end
            if ~isfield(self.complementary_information, 'long_run_file')
                self.complementary_information.long_run_file = '—';
            end

        end


        function add_var_tab_2_inputs(self)
            % initiate lines
            lines = string([]);
            % header for tab 2
            lines = [lines;'Specification'];
            lines = [lines;'-----------------'];
            lines = [lines;' '];
            % VAR type
            var_type = self.complementary_information.model_type;
            if var_type == 1
                model = 'Maximum Likelihood VAR';
            elseif var_type == 2
                model = 'Minnesota Bayesian VAR';          
            elseif var_type == 3
                model = 'Normal-Wishart Bayesian VAR';        
            elseif var_type == 4
                model = 'Independent Bayesian VAR';
            elseif var_type == 5
                model = 'Dummy Observation Bayesian VAR';             
            elseif var_type == 6
                model = 'Large Bayesian VAR';
            elseif var_type == 7
                model = 'Bayesian Proxy-SVAR';   
            end
            lines = [lines;'VAR type: ' model];
            % burn-in
            if ismember(var_type, [4 6 7])
                burnin = num2str(self.model.burnin);
                lines = [lines;'burn-in: ' burnin];
            end
            % iterations
            if var_type ~= 1
                iterations = num2str(self.model.iterations);
                lines = [lines;'iterations: ' iterations];
            end
            % credibility level
            model_credibility = num2str(self.model.credibility_level);
            lines = [lines;'credibility level: ' model_credibility];
            % constant, trend and quadratic trend
            constant = cu.bool_to_string(self.model.constant);
            lines = [lines;'constant: ' constant];      
            trend = cu.bool_to_string(self.model.trend);
            lines = [lines;'trend: ' trend]; 
            quadratic_trend = cu.bool_to_string(self.model.quadratic_trend);
            lines = [lines;'quadratic trend: ' quadratic_trend];
            % hyperparameters
            lags = num2str(self.model.p);
            lines = [lines;'lags: ' lags];
            if var_type ~= 1
                if numel(self.model.ar_coefficients) == 1
                    ar_coefficients = num2str(self.model.ar_coefficients);
                else
                    ar_coefficients = iu.array_to_char(self.model.ar_coefficients);
                end
                lines = [lines;'AR coefficients: ' ar_coefficients];
                pi1 = num2str(self.model.pi1);
                lines = [lines;'pi1 (overall tightness): ' pi1];
                if ismember(var_type, [2 4 6])
                    pi2 = num2str(self.model.pi2);
                    lines = [lines;'pi2 (cross-variable shrinkage): ' pi2];
                end
                pi3 = num2str(self.model.pi3);
                lines = [lines;'pi3 (lag decay): ' pi3];
                pi4 = num2str(self.model.pi4);
                lines = [lines;'pi4 (exogenous slackness): ' pi4];
                if var_type ~= 7
                    pi5 = num2str(self.model.pi5);
                    lines = [lines;'pi5 (sums-of-coefficients tightness): ' pi5];
                    pi6 = num2str(self.model.pi6);
                    lines = [lines;'pi6 (initial observation tightness): ' pi6];
                    pi7 = num2str(self.model.pi7);
                    lines = [lines;'pi7 (long-run tightness): ' pi7];
                end
            end
            if var_type == 7
                proxy_variables = iu.array_to_char(self.complementary_information.proxy_variables);
                lines = [lines;'proxy variables: ' proxy_variables];
                lamda = num2str(self.model.lamda);
                lines = [lines;'lambda (relevance): ' lamda];
                proxy_prior = self.model.proxy_prior;
                if proxy_prior == 1
                    prior_scheme = 'uninformative';
                elseif proxy_prior == 2
                    prior_scheme = 'Minnesota';
                end
                lines = [lines;'prior scheme: ' prior_scheme];
            end
            % VAR options: in-sample fit
            if islogical(self.complementary_information.insample_fit')
                insample_fit = cu.bool_to_string(self.complementary_information.insample_fit); 
            else
                insample_fit = self.complementary_information.insample_fit;
            end
            lines = [lines;'in-sample fit: ' insample_fit];
            % VAR options: constrained coefficients
            if ismember(var_type, [2 4 6])
                constrained_coefficients = cu.bool_to_string(self.model.constrained_coefficients); 
                lines = [lines;'constrained coefficients: ' constrained_coefficients];
            end
            % VAR options: sums-of-coefficients
            if ismember(var_type, [2 3 4 5 6])
                sums_of_coefficients = cu.bool_to_string(self.model.sums_of_coefficients);   
                lines = [lines;'sums-of-coefficients: ' sums_of_coefficients];
            end
            % VAR options: dummy initial observation
            if ismember(var_type, [2 3 4 5 6])
                initial_observation = cu.bool_to_string(self.model.dummy_initial_observation);
                lines = [lines;'dummy initial observation: ' initial_observation];
            end
            % VAR options: long-run prior
            if ismember(var_type, [2 3 4 5 6])
                long_run_prior = cu.bool_to_string(self.model.long_run_prior);  
                lines = [lines;'long-run prior: ' long_run_prior];
            end
            % VAR options: stationary prior
            if ismember(var_type, [2 3 4 5 6])
                stationary_prior = cu.bool_to_string(self.model.stationary_prior); 
                lines = [lines;'stationary prior: ' stationary_prior];
            end
            % VAR options: marginal likelihood
            if ismember(var_type, [2 3 4])
                if islogical(self.complementary_information.marginal_likelihood)
                    marginal_likelihood = cu.bool_to_string(self.complementary_information.marginal_likelihood);   
                else
                    marginal_likelihood = self.complementary_information.marginal_likelihood;
                end
                lines = [lines;'marginal likelihood: ' marginal_likelihood];
            end
            % VAR options: hyperparameter optimization
            if var_type == 2 || var_type == 3
                hyperparameter_optimization = cu.bool_to_string(self.model.hyperparameter_optimization);
                lines = [lines;'hyperparameter optimization: ' hyperparameter_optimization];
            end
            % VAR options: constrained coefficient file
            if ismember(var_type, [2 4 6])
                constrained_coefficient_file = self.complementary_information.coefficients_file;
                lines = [lines;'constrained coefficients file: ' constrained_coefficient_file];
            end
            % VAR options: long-run prior file
            if ismember(var_type, [2 3 4 5 6])
                long_run_file = self.complementary_information.long_run_file;
                lines = [lines;'long-run prior file: ' long_run_file];
            end
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end


        function add_var_tab_3_inputs(self)
            % initiate lines
            lines = string([]);
            % header for tab 1
            lines = [lines;'Applications'];
            lines = [lines;'---------'];
            lines = [lines;' '];
            % forecasts
            if islogical(self.complementary_information.forecast)
                forecast = cu.bool_to_string(self.complementary_information.forecast);
            else
                forecast = self.complementary_information.forecast;
            end
            lines = [lines;'forecast: ' forecast];
            forecast_credibility = num2str(self.complementary_information.forecast_credibility);
            lines = [lines;'credibility level, forecasts: ' forecast_credibility];
            % conditional forecasts
            if islogical(self.complementary_information.conditional_forecast)
                conditional_forecast = cu.bool_to_string(self.complementary_information.conditional_forecast);
            else
                conditional_forecast = self.complementary_information.conditional_forecast;
            end
            lines = [lines;'conditional forecast: ' conditional_forecast];
            conditional_forecast_credibility = num2str(self.complementary_information.conditional_forecast_credibility);
            lines = [lines;'credibility level, conditional forecasts: ' conditional_forecast_credibility];
            % impulse response function
            if islogical(self.complementary_information.irf)
                irf = cu.bool_to_string(self.complementary_information.irf);
            else
                irf = self.complementary_information.irf;
            end
            lines = [lines;'impulse response function: ' irf]; 
            irf_credibility = num2str(self.complementary_information.irf_credibility);
            lines = [lines;'credibility level, impulse response function: ' irf_credibility];
            % forecast error variance decomposition
            if islogical(self.complementary_information.fevd)
                fevd = cu.bool_to_string(self.complementary_information.fevd);
            else
                fevd = self.complementary_information.fevd;
            end
            lines = [lines;'forecast error variance decomposition: ' fevd];
            fevd_credibility = num2str(self.complementary_information.fevd_credibility);
            lines = [lines;'credibility level, forecast error variance decomposition: ' fevd_credibility];
            % historical decomposition
            if islogical(self.complementary_information.hd)
                hd = cu.bool_to_string(self.complementary_information.hd);
            else
                hd = self.complementary_information.hd;
            end
            lines = [lines;'historical decomposition: ' hd];
            hd_credibility = num2str(self.complementary_information.hd_credibility);
            lines = [lines;'credibility level, historical decomposition: ' hd_credibility];
            % forecast periods
            if isnumeric(self.complementary_information.forecast_periods)
                forecast_periods = num2str(self.complementary_information.forecast_periods);
            else
                forecast_periods = self.complementary_information.forecast_periods;
            end
            lines = [lines;'forecast periods: ' forecast_periods];
            % conditional forecast type
            if isnumeric(self.complementary_information.conditional_forecast_type)
                if self.complementary_information.conditional_forecast_type == 1
                    conditional_forecast_type = 'agnostic';
                elseif self.complementary_information.conditional_forecast_type == 2
                    conditional_forecast_type = 'structural shocks';
                end
            else
                conditional_forecast_type = self.complementary_information.conditional_forecast_type;
            end
            lines = [lines;'conditional forecast type: ' conditional_forecast_type];
            % forecast file
            forecast_file = self.complementary_information.forecast_file;
            lines = [lines;'forecast file: ' forecast_file];
            % conditional forecast file
            conditional_forecast_file = self.complementary_information.conditional_forecast_file;
            lines = [lines;'conditional forecast file: ' conditional_forecast_file];
            % forecast evaluation
            if islogical(self.complementary_information.forecast_evaluation)
                forecast_evaluation = cu.bool_to_string(self.complementary_information.forecast_evaluation);
            else
                forecast_evaluation = self.complementary_information.forecast_evaluation;
            end
            lines = [lines;'forecast evaluation: ' forecast_evaluation];
            % irf periods
            if isnumeric(self.complementary_information.irf_periods)
                irf_periods = num2str(self.complementary_information.irf_periods);
            else
                irf_periods = self.complementary_information.irf_periods;
            end
            lines = [lines;'IRF periods: ' irf_periods];
            % structural identification
            if isnumeric(self.complementary_information.structural_identification)
                structural_identification = num2str(self.complementary_information.structural_identification);
            else
                structural_identification = self.complementary_information.structural_identification;
            end
            lines = [lines;'structural identification: ' structural_identification];
            % structural identification file        
            structural_identification_file = self.complementary_information.structural_identification_file;
            lines = [lines;'structural identification file: ' structural_identification_file];
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end


        function make_var_summary(self)
            % initiate string list
            self.estimation_summary = string([]);
            % add model header
            self.add_var_header();
            % add estimation header
            self.add_var_estimation_header();
            % make list of regressors
            self.regressors_and_index();
            % loop over equations
            for i=1:self.model.n
                % add coefficient summary
                self.add_var_coefficient_summary(i);
                % residual and shock variance
                self.add_shock_variance_summary(i);
                % in-sample fit criteria
                self.add_var_insample_evaluation(i);
            end
            % residual variance-covariance matrix
            self.add_residual_matrix_summary();
            % structural shocks variance-covariance matrix
            self.add_shock_matrix_summary();
            % structural identification matrix
            self.add_structural_identification_matrix_summary();
            % add forecast evaluation criteria, if relevant
            self.add_var_forecast_evaluation();
        end
        

        function make_var_application_summary(self)
            % in-sample fit measures
            self.make_var_insample_fit_summary();
            % forecasts
            self.make_var_forecast_summary();        
            % conditional forecasts
            self.make_var_conditional_forecast_summary();  
            % impulse response function
            self.make_var_irf_summary();        
            % forecast error variance decomposition
            self.make_var_fevd_summary();    
            % historical decomposition
            self.make_var_hd_summary();
        end


        function save_var_application(self, path)
            % save in-sample fit
            self.save_var_insample_fit_summary(path);
            % save forecasts
            self.save_var_forecast_summary(path);   
            % save conditional forecasts
            self.save_var_conditional_forecast_summary(path);      
            % save impulse response function
            self.save_var_irf_summary(path);
            % save forecast error variance decomposition
            self.save_var_fevd_summary(path);
            % save historical decomposition
            self.save_var_hd_summary(path);
        end


        function add_var_header(self)
            % recover model name and create header
            model_name = self.complementary_information.model_name;
            self.estimation_summary = [self.estimation_summary;cu.model_header(model_name)];
        end


        function add_var_estimation_header(self)
            % initiate lines
            lines = string([]);
            % first row: estimation sample and estimation start
            sample_start = self.complementary_information.sample_start;
            sample_end = self.complementary_information.sample_end;
            if isempty(sample_start) || isempty(sample_end)
                sample = '—';
            else
                sample = [sample_start '  ' sample_end];
            end
            estimation_start = self.complementary_information.estimation_start;
            left_element = sprintf(['Sample:' sprintf('%+31s', sample)]);
            right_element = sprintf(['Est. start:' sprintf('%+27s', estimation_start)]);
            lines = [lines; [left_element '    ' right_element]];
            % second row: observations and estimation complete   
            T = num2str(self.model.T);
            estimation_end = self.complementary_information.estimation_end;
            left_element = sprintf(['No. observations:' sprintf('%+21s', T)]);
            right_element = sprintf(['Est. complete:' sprintf('%+24s', estimation_end)]);
            lines = [lines; [left_element '    ' right_element]];
            % third row: frequency and lags
            frequency = self.complementary_information.frequency;
            lags = num2str(self.model.p);
            left_element = sprintf(['Frequency:' sprintf('%+28s', frequency)]);
            right_element = sprintf(['Lags:' sprintf('%+33s', lags)]);
            lines = [lines; [left_element '    ' right_element]];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function regressors_and_index(self)
            endogenous = self.complementary_information.endogenous_variables;
            exogenous = self.complementary_information.exogenous_variables;
            constant = self.model.constant;
            trend = self.model.trend;
            quadratic_trend = self.model.quadratic_trend;
            n = self.model.n;
            m = self.model.m;
            p = self.model.p;
            k = self.model.k;
            regressors = cu.make_regressors(endogenous, exogenous, constant, trend, quadratic_trend, n, p);
            coefficient_index = cu.make_index(n, m, p, k);
            self.regressors = regressors;
            self.coefficient_index = coefficient_index;
        end

        
        function add_var_coefficient_summary(self, i)
            lines = string([]);
            endogenous_variables = self.complementary_information.endogenous_variables;
            credibility_level = self.model.credibility_level;
            lines = [lines;cu.equation_header(['Equation: ' char(endogenous_variables(i))])];
            lines = [lines;cu.coefficient_header(credibility_level)];
            lines = [lines;cu.string_line('VAR coefficients beta:')];
            % loop over equation coefficients
            coefficient_index = self.coefficient_index;
            regressors = self.regressors;
            for j=1:self.model.k
                regressor = char(regressors(j));
                index = coefficient_index(j);
                coefficient = self.model.beta_estimates(index,i,1);
                standard_deviation = self.model.beta_estimates(index,i,4);
                lower_bound = self.model.beta_estimates(index,i,2);
                upper_bound = self.model.beta_estimates(index,i,3);
                lines = [lines;cu.string_line(cu.parameter_estimate_line(regressor,...
                         coefficient, standard_deviation, lower_bound, upper_bound))];
            end
            lines = [lines;cu.hyphen_dashed_line()];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_shock_variance_summary(self, i)
            lines = string([]);
            residual_variance = self.model.Sigma_estimates(i,i);
            if ~isempty(self.model.Gamma_estimates)
                shock_variance = self.model.Gamma_estimates(i);
            else
                shock_variance = '';
            end
            lines = [lines;cu.variance_line(residual_variance, shock_variance)];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_var_insample_evaluation(self, i)
            % initiate lines
            lines = string([]);
            % check if in-sample evaluation has been conducted
            if ~isempty(self.model.insample_evaluation)
                lines = [lines;cu.hyphen_dashed_line()];
                ssr = self.model.insample_evaluation.ssr(i);
                r2 = self.model.insample_evaluation.r2(i);
                adj_r2 = self.model.insample_evaluation.adj_r2(i); 
                model_type = self.complementary_information.model_type;
                if model_type == 1
                    aic = self.model.insample_evaluation.aic;
                    bic = self.model.insample_evaluation.bic;
                else
                    aic = [];
                    bic = [];
                end
                if isprop(self.model,'m_y') && ~isempty(self.model.m_y)
                    m_y = self.model.m_y;
                else
                    m_y = [];
                end
                lines = [lines;cu.insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic)];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function add_residual_matrix_summary(self)
            Sigma = self.model.Sigma_estimates;
            n = self.model.n;
            endogenous_variables = self.complementary_information.endogenous_variables;  
            lines = string([]);
            lines = [lines;cu.equation_header('Residual variance-covariance Sigma')];
            lines = [lines;cu.variance_covariance_summary(Sigma, n, endogenous_variables, 'var.Sigma_estimates')];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_shock_matrix_summary(self)
            if isprop(self.model, 'Gamma_estimates') && ~isempty(self.model.Gamma_estimates)
                Gamma = diag(self.model.Gamma_estimates);
                n = self.model.n;
                endogenous_variables = self.complementary_information.endogenous_variables; 
                lines = string([]);
                lines = [lines;cu.intermediate_header('Structural shocks variance-covariance Gamma')];
                lines = [lines;cu.variance_covariance_summary(Gamma, n, endogenous_variables, 'var.Gamma_estimates')];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function add_structural_identification_matrix_summary(self)
            lines = string([]);
            if isprop(self.model, 'H_estimates') && ~isempty(self.model.H_estimates)
                H = self.model.H_estimates;
                n = self.model.n;
                endogenous_variables = self.complementary_information.endogenous_variables;       
                lines = [lines;cu.intermediate_header('Structural identification matrix H')];
                lines = [lines;cu.variance_covariance_summary(H, n, endogenous_variables, 'var.H_estimates')];
            end
            self.estimation_summary = [self.estimation_summary;lines]; 
        end


        function add_var_forecast_evaluation(self)
            lines = string([]);
            if ~isempty(self.model.forecast_evaluation_criteria) 
                endogenous_variables = self.complementary_information.endogenous_variables;  
                forecast_evaluation_criteria = self.model.forecast_evaluation_criteria;
                % regular forecast evaluation criteria
                lines = [lines;cu.equation_header('Forecast evaluation criteria')];
                lines = [lines;'                 RMSE        MAE       MAPE    Theil-U       Bias               '];
                rmse = forecast_evaluation_criteria.rmse;
                mae = forecast_evaluation_criteria.mae;
                mape = forecast_evaluation_criteria.mape;
                theil_u = forecast_evaluation_criteria.theil_u;
                bias = forecast_evaluation_criteria.bias;
                for i=1:self.model.n
                    lines = [lines;cu.forecast_evaluation_line(endogenous_variables(i), ...
                                 rmse(i), mae(i), mape(i), theil_u(i), bias(i))];
                end
                % Bayesian criteria: log score
                if isfield(forecast_evaluation_criteria, 'log_score')
                    lines = [lines;cu.intermediate_header('Log score')];
                    log_score = forecast_evaluation_criteria.log_score;
                    joint_log_score = forecast_evaluation_criteria.joint_log_score;
                    lines = [lines;cu.forecast_evaluation_summary(log_score, joint_log_score, ...
                             endogenous_variables, 'var.forecast_evaluation_criteria[''log_score'']')];
                end
                % Bayesian criteria: CRPS
                if isfield(forecast_evaluation_criteria, 'crps')
                    lines = [lines;cu.intermediate_header('CRPS')];
                    crps = forecast_evaluation_criteria.crps;
                    joint_crps = forecast_evaluation_criteria.joint_crps;
                    lines = [lines;cu.forecast_evaluation_summary(crps, joint_crps, ...
                             endogenous_variables, 'var.forecast_evaluation_criteria[''crps'']')];
                end
            end
            lines = [lines;cu.equal_dashed_line()];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function make_var_insample_fit_summary(self)
            % run only if in-sample fit has been run
            if ~isempty(self.model.fitted_estimates)
                Y = self.model.Y;
                endogenous_variables = self.complementary_information.endogenous_variables;
                fitted = self.model.fitted_estimates;
                residuals = self.model.residual_estimates;
                n = self.model.n;
                p = self.model.p;
                index = self.complementary_information.dates(p+1:end);
                fitted_table = table(index);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    header = [string([variable '_actual']) [variable '_fit_med'] [variable '_fit_low'] ...
                             [variable '_fit_upp'] [variable '_res_med'] [variable '_res_low'] [variable '_res_upp']];
                    actual = Y(:,i);
                    fit_med = fitted(:,i,1);
                    fit_low = fitted(:,i,2);
                    fit_upp = fitted(:,i,3);
                    res_med = residuals(:,i,1);
                    res_low = residuals(:,i,2);
                    res_upp = residuals(:,i,3);   
                    variable_table = table(actual, fit_med, fit_low, fit_upp, res_med, res_low, res_upp);
                    variable_table.Properties.VariableNames = header;
                    fitted_table = [fitted_table variable_table];
                end
                self.application_summary.insample_fit = fitted_table;
            end
        end


        function make_var_forecast_summary(self)
            % run only if forecast has been run
            if ~isempty(self.model.forecast_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                p = self.model.p;
                Y = self.model.Y;
                insample_index = self.complementary_information.dates(p+1:end);
                forecasts = self.model.forecast_estimates;
                forecast_index = self.complementary_information.forecast_dates;
                forecast_table = table([insample_index;forecast_index]);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    header = [string([variable '_actual']) [variable '_med'] [variable '_low'] [variable '_fit_upp']];
                    insample_table = table(Y(:,i));
                    insample_table.Var2(:) = NaN;
                    insample_table.Var3(:) = NaN;
                    insample_table.Var4(:) = NaN;
                    insample_table(end,:) = insample_table(end,1);
                    Var2  = forecasts(:,i,1);
                    Var3  = forecasts(:,i,2);
                    Var4  = forecasts(:,i,3);
                    prediction_table = table(Var2, Var3, Var4);
                    prediction_table.Var1(:) = NaN;
                    variable_table = [insample_table;prediction_table];
                    variable_table.Properties.VariableNames = header;
                    forecast_table = [forecast_table variable_table];
                end
                self.application_summary.forecast = forecast_table;
            end
        end


        function make_var_conditional_forecast_summary(self)
            % run only if conditional forecast has been run
            if isprop(self.model, 'conditional_forecast_estimates') && ~isempty(self.model.conditional_forecast_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                p = self.model.p;
                Y = self.model.Y;
                insample_index = self.complementary_information.dates(p+1:end);
                if ~isempty(self.model.conditional_forecast_estimates)
                    forecasts = self.model.conditional_forecast_estimates;
                end
                forecast_index = self.complementary_information.forecast_dates;
                forecast_table = table([insample_index;forecast_index]);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    header = [string([variable '_actual']) [variable '_med'] [variable '_low'] [variable '_fit_upp']];
                    insample_table = table(Y(:,i));
                    insample_table.Var2(:) = NaN;
                    insample_table.Var3(:) = NaN;
                    insample_table.Var4(:) = NaN;
                    insample_table(end,:) = insample_table(end,1);
                    Var2  = forecasts(:,i,1);
                    Var3  = forecasts(:,i,2);
                    Var4  = forecasts(:,i,3);
                    prediction_table = table(Var2, Var3, Var4);
                    prediction_table.Var1(:) = NaN;
                    variable_table = [insample_table;prediction_table];
                    variable_table.Properties.VariableNames = header;
                    forecast_table = [forecast_table variable_table];
                end
                self.application_summary.conditional_forecast = forecast_table;
            end
        end


        function make_var_irf_summary(self)
            endogenous_variables = self.complementary_information.endogenous_variables;
            n = self.model.n;
            % run only if IRF has been run
            if ~isempty(self.model.irf_estimates)
                irf = permute(self.model.irf_estimates, [3 4 1 2]);
                index = (1:size(irf,1))';
                irf_table = table(index);
                for i=1:n
                    for j=1:n
                        variable = char(endogenous_variables(i));
                        shock = ['shock' num2str(j)];
                        header = [string([variable '_' shock '_med']) [variable '_' shock '_low'] ...
                                 [variable '_' shock '_upp']];
                        Var1 = irf(:,1,i,j);
                        Var2 = irf(:,2,i,j);
                        Var3 = irf(:,3,i,j);
                        variable_table = table(Var1, Var2, Var3);
                        variable_table.Properties.VariableNames = header;
                        irf_table = [irf_table variable_table];
                    end
                end
                self.application_summary.irf = irf_table;
            end
            % run only if exogenous IRF have been computed
            if isprop(self.model, 'exo_irf_estimates') && ~isempty(self.model.exo_irf_estimates)
                exo_irf = permute(self.model.exo_irf_estimates, [3 4 1 2]);
                index = (1:size(irf,1))';
                exo_irf_table = table(index);
                exogenous_variables = self.complementary_information.exogenous_variables;
                n_exo = numel(exogenous_variables);
                for i=1:n
                    for j=1:n_exo
                        variable = char(endogenous_variables(i));
                        shock = char(exogenous_variables(j));
                        header = [string([variable '_' shock '_med']) [variable '_' shock '_low'] ...
                                 [variable '_' shock '_upp']];
                        Var1 = exo_irf(:,1,i,j);
                        Var2 = exo_irf(:,2,i,j);
                        Var3 = exo_irf(:,3,i,j);
                        variable_table = table(Var1, Var2, Var3);
                        variable_table.Properties.VariableNames = header;
                        exo_irf_table = [exo_irf_table variable_table];
                    end
                end
                self.application_summary.exo_irf = exo_irf_table;
            end
        end


        function make_var_fevd_summary(self)
            % run only if FEVD has been run
            if isprop(self.model, 'fevd_estimates') && ~isempty(self.model.fevd_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                fevd = permute(self.model.fevd_estimates, [3 4 1 2]);
                index = (1:size(fevd,1))';
                fevd_table = table(index);
                for i=1:n
                    for j=1:n
                        variable = char(endogenous_variables(i));
                        shock = ['shock' num2str(j)];
                        header = [string([variable '_' shock '_med']) [variable '_' shock '_low'] ...
                                 [variable '_' shock '_upp']];
                        Var1 = fevd(:,1,i,j);
                        Var2 = fevd(:,2,i,j);
                        Var3 = fevd(:,3,i,j);
                        variable_table = table(Var1, Var2, Var3);
                        variable_table.Properties.VariableNames = header;
                        fevd_table = [fevd_table variable_table];
                    end
                end
                self.application_summary.fevd = fevd_table;
            end
        end


        function make_var_hd_summary(self)
            % run only if HD has been run
            if isprop(self.model, 'hd_estimates') && ~isempty(self.model.hd_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                p = self.model.p;
                hd = permute(self.model.hd_estimates, [3 4 1 2]);
                index = self.complementary_information.dates(p+1:end);
                hd_table = table(index);
                for i=1:n
                    for j=1:n
                        variable = char(endogenous_variables(i));
                        shock = ['shock' num2str(j)];
                        header = [string([variable '_' shock '_med']) [variable '_' shock '_low'] ...
                                 [variable '_' shock '_upp']];
                        Var1 = hd(:,1,i,j);
                        Var2 = hd(:,2,i,j);
                        Var3 = hd(:,3,i,j);
                        variable_table = table(Var1, Var2, Var3);
                        variable_table.Properties.VariableNames = header;
                        hd_table = [hd_table variable_table];
                    end
                end
                self.application_summary.hd = hd_table;
            end
        end


        function save_var_insample_fit_summary(self, path)
            if isfield(self.application_summary,'insample_fit')
                insample_fit_summary = self.application_summary.insample_fit;
                full_path = fullfile(path, 'insample_fit.csv');
                writetable(insample_fit_summary, full_path);
            end
        end


        function save_var_forecast_summary(self, path)
            if isfield(self.application_summary,'forecast')
                forecast_summary = self.application_summary.forecast;
                full_path = fullfile(path, 'forecast.csv');
                writetable(forecast_summary, full_path);
            end
        end


        function save_var_conditional_forecast_summary(self, path)
            if isfield(self.application_summary,'conditional_forecast')
                conditional_forecast_summary = self.application_summary.conditional_forecast;
                full_path = fullfile(path, 'conditional_forecast.csv');
                writetable(conditional_forecast_summary, full_path);
            end
        end


        function save_var_irf_summary(self, path)
            if isfield(self.application_summary,'irf')
                irf_summary = self.application_summary.irf;
                full_path = fullfile(path, 'irf.csv');
                writetable(irf_summary, full_path);
            end
        end


        function save_var_fevd_summary(self, path)
            if isfield(self.application_summary,'fevd')
                fevd_summary = self.application_summary.fevd;
                full_path = fullfile(path, 'fevd.csv');
                writetable(fevd_summary, full_path);
            end
        end


        function save_var_hd_summary(self, path)
            if isfield(self.application_summary,'hd')
                hd_summary = self.application_summary.hd;
                full_path = fullfile(path, 'hd.csv');
                writetable(hd_summary, full_path);
            end
        end

    end

end
