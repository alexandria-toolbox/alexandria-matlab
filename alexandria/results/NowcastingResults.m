classdef NowcastingResults < handle
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        mfbvar_regressors
        mfbvar_coefficient_index
        bdfm_loadings_regressors
        bdfm_factor_regressors
        bdfm_residual_regressors
        bdfm_loadings_index
        bdfm_factor_index
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    


        function self = NowcastingResults()          
        end


    end
    
    
    methods (Access = protected, Hidden = true)
        

        function complete_nowcasting_information(self)
            % endogenous and exogenous variables
            if ~isfield(self.complementary_information, 'endogenous_variables')
                if self.complementary_information.model_type == 3
                    self.complementary_information.endogenous_variables = "y1";
                else
                    n_endo = size(self.model.endogenous,2);
                    self.complementary_information.endogenous_variables = "y"+(1:n_endo);
                end
            end
            if ~isfield(self.complementary_information, 'exogenous_variables')
                if self.complementary_information.model_type == 2 || isempty(self.model.exogenous)
                    self.complementary_information.exogenous_variables = "none";
                else
                    n_exo = size(self.model.exogenous,2);
                    self.complementary_information.exogenous_variables = "x"+(1:n_exo);
                end
            end
            % sample dates
            if ~isfield(self.complementary_information, 'dates')
                T = self.model.T;
                if self.complementary_information.model_type == 2
                    p = 0;
                else
                    p = self.model.p;
                end
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
                if isprop(self.model, 'conditional_forecast_estimates') && ~isempty(self.model.conditional_forecast_estimates)
                    f_periods = size(self.model.conditional_forecast_estimates,1);
                    T = self.model.T;
                    self.complementary_information.conditional_forecast_dates = (T+1:T+f_periods)';
                else
                    self.complementary_information.conditional_forecast_dates = [];
                end
            end
        end


        function add_nowcasting_tab_2_inputs(self)
            % initiate lines
            lines = string([]);
            % header for tab 2
            lines = [lines;'Specification'];
            lines = [lines;'-----------------'];
            lines = [lines;' '];
            % model
            model = self.complementary_information.model_name;  
            lines = [lines;'model: ' model];
            % iterations
            iterations = num2str(self.model.iterations);
            lines = [lines;'iterations: ' iterations];
            % burn-in
            burnin = num2str(self.model.burnin);
            lines = [lines;'burn-in: ' burnin];
            % credibility level
            model_credibility = num2str(self.model.credibility_level);
            lines = [lines;'credibility level: ' model_credibility];
            % get model type, midas regression, mfbvar or bdfm
            model_type = self.complementary_information.model_type;        
            % mfbvar hyperparmeters
            if model_type == 1
                % constant, trend and quadratic trend
                constant = cu.bool_to_string(self.model.constant);
                lines = [lines;'constant: ' constant];      
                trend = cu.bool_to_string(self.model.trend);
                lines = [lines;'trend: ' trend]; 
                quadratic_trend = cu.bool_to_string(self.model.quadratic_trend);
                lines = [lines;'quadratic trend: ' quadratic_trend];
                % decomposition
                decomposition = cu.bool_to_string(self.model.decomposition);
                lines = [lines;'decomposition: ' decomposition]; 
                % hyperparameters
                lags = num2str(self.model.p);
                lines = [lines;'lags: ' lags];                
                if numel(self.model.ar_coefficients) == 1
                    ar_coefficients = num2str(self.model.ar_coefficients);
                else
                    ar_coefficients = iu.array_to_char(self.model.ar_coefficients);
                end
                lines = [lines;'AR coefficients: ' ar_coefficients];
                pi1 = num2str(self.model.pi1);
                lines = [lines;'pi1 (overall tightness): ' pi1];
                pi2 = num2str(self.model.pi2);
                lines = [lines;'pi2 (cross-variable shrinkage): ' pi2];
                pi3 = num2str(self.model.pi3);
                lines = [lines;'pi3 (lag decay): ' pi3];
                pi4 = num2str(self.model.pi4);
                lines = [lines;'pi4 (exogenous slackness): ' pi4];
                decomposition_file = self.complementary_information.decomposition_file;
                lines = [lines;'decomposition file: ' decomposition_file];
            % bdfm hyperparameters
            elseif model_type == 2
                m = num2str(self.model.m);
                lines = [lines;'m (factors): ' m];                
                q = num2str(self.model.q);
                lines = [lines;'q (loadings lags): ' q]; 
                p = num2str(self.model.p);
                lines = [lines;'p (factor lags): ' p]; 
                r = num2str(self.model.r);
                lines = [lines;'r (residual lags): ' r]; 
                sigma = num2str(self.model.sigma);
                lines = [lines;'sigma (residual variance): ' sigma];
                omega = num2str(self.model.omega);
                lines = [lines;'omega (factor variance): ' omega];
                delta1 = num2str(self.model.delta1);
                lines = [lines;'delta1 (loadings tightness): ' delta1];
                pi1 = num2str(self.model.pi1);
                lines = [lines;'pi1 (factor tightness): ' pi1];
                pi2 = num2str(self.model.pi2);
                lines = [lines;'pi2 (cross-variable shrinkage): ' pi2];
                pi3 = num2str(self.model.pi3);
                lines = [lines;'pi3 (lag decay): ' pi3];
                omega1 = num2str(self.model.omega1);
                lines = [lines;'omega1 (residual tightness): ' omega1];
            % midas regression hyperparameters
            elseif model_type == 3
                representation = self.model.representation;
                lines = [lines;'representation: ' representation];                
                prior_type = self.model.prior_type;
                lines = [lines;'prior: ' prior_type];
                endogenous_lags = num2str(self.model.endogenous_lags);
                lines = [lines;'endogenous lags: ' endogenous_lags];   
                if numel(self.model.exogenous_lags) == 1
                    exogenous_lags = num2str(self.model.exogenous_lags);
                else
                    exogenous_lags = iu.array_to_char(self.model.exogenous_lags);
                end
                lines = [lines;'exogenous lags: ' exogenous_lags];  
                polynomial_order = num2str(self.model.polynomial_order);
                lines = [lines;'polynomial order: ' polynomial_order];  
                omega1 = num2str(self.model.omega1);
                lines = [lines;'omega1 (endogenous tightness): ' omega1];
                omega2 = num2str(self.model.omega2);
                lines = [lines;'omega2 (endogenous lag decay): ' omega2];
                upsilon1 = num2str(self.model.upsilon1);
                lines = [lines;'upsilon1 (exogenous tightness): ' upsilon1];
                upsilon2 = num2str(self.model.upsilon2);
                lines = [lines;'upsilon2 (exogenous lag decay): ' upsilon2];
            end
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end


        function add_nowcasting_tab_3_inputs(self)
            % get model type, midas regression, mfbvar or bdfm
            model_type = self.complementary_information.model_type;   
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
            if model_type == 1
                if islogical(self.complementary_information.conditional_forecast)
                    conditional_forecast = cu.bool_to_string(self.complementary_information.conditional_forecast);
                else
                    conditional_forecast = self.complementary_information.conditional_forecast;
                end
                lines = [lines;'conditional forecast: ' conditional_forecast];
                conditional_forecast_credibility = num2str(self.complementary_information.conditional_forecast_credibility);
                lines = [lines;'credibility level, conditional forecasts: ' conditional_forecast_credibility];
            end
            % impulse response function
            if model_type == 1 || model_type == 2
                if islogical(self.complementary_information.irf)
                    irf = cu.bool_to_string(self.complementary_information.irf);
                else
                    irf = self.complementary_information.irf;
                end
                lines = [lines;'impulse response function: ' irf]; 
                irf_credibility = num2str(self.complementary_information.irf_credibility);
                lines = [lines;'credibility level, impulse response function: ' irf_credibility];
            end
            % forecast error variance decomposition
            if model_type == 1 || model_type == 2
                if islogical(self.complementary_information.fevd)
                    fevd = cu.bool_to_string(self.complementary_information.fevd);
                else
                    fevd = self.complementary_information.fevd;
                end
                lines = [lines;'forecast error variance decomposition: ' fevd];
                fevd_credibility = num2str(self.complementary_information.fevd_credibility);
                lines = [lines;'credibility level, forecast error variance decomposition: ' fevd_credibility];
            end
            % historical decomposition
            if model_type == 1 || model_type == 2
                if islogical(self.complementary_information.hd)
                    hd = cu.bool_to_string(self.complementary_information.hd);
                else
                    hd = self.complementary_information.hd;
                end
                lines = [lines;'historical decomposition: ' hd];
                hd_credibility = num2str(self.complementary_information.hd_credibility);
                lines = [lines;'credibility level, historical decomposition: ' hd_credibility];
            end
            % forecast periods
            if isnumeric(self.complementary_information.forecast_periods)
                forecast_periods = num2str(self.complementary_information.forecast_periods);
            else
                forecast_periods = self.complementary_information.forecast_periods;
            end
            lines = [lines;'forecast periods: ' forecast_periods];
            % conditional forecast type
            if model_type == 1
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
            end
            % forecast file
            forecast_file = self.complementary_information.forecast_file;
            lines = [lines;'forecast file: ' forecast_file];
            % conditional forecast file
            if model_type == 1
                conditional_forecast_file = self.complementary_information.conditional_forecast_file;
                lines = [lines;'conditional forecast file: ' conditional_forecast_file];
            end
            % forecast evaluation
            if islogical(self.complementary_information.forecast_evaluation)
                forecast_evaluation = cu.bool_to_string(self.complementary_information.forecast_evaluation);
            else
                forecast_evaluation = self.complementary_information.forecast_evaluation;
            end
            lines = [lines;'forecast evaluation: ' forecast_evaluation];
            % irf periods
            if model_type == 1 || model_type == 2
                if isnumeric(self.complementary_information.irf_periods)
                    irf_periods = num2str(self.complementary_information.irf_periods);
                else
                    irf_periods = self.complementary_information.irf_periods;
                end
                lines = [lines;'IRF periods: ' irf_periods];
            end
            % structural identification
            if model_type == 1
                if isnumeric(self.complementary_information.structural_identification)
                    structural_identification = num2str(self.complementary_information.structural_identification);
                else
                    structural_identification = self.complementary_information.structural_identification;
                end
                lines = [lines;'structural identification: ' structural_identification];
            end
            % structural identification file  
            if model_type == 1
                structural_identification_file = self.complementary_information.structural_identification_file;
                lines = [lines;'structural identification file: ' structural_identification_file];
            end
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end


        function make_nowcasting_summary(self)
            % initiate string list
            self.estimation_summary = string([]);
            % add model header
            self.add_nowcasting_header();
            % add estimation header
            self.add_nowcasting_estimation_header();
            % get model type, midas regression, mfbvar or bdfm
            model_type = self.complementary_information.model_type;   
            % mfbvar summary
            if model_type == 1
                self.add_mfbvar_summary();
            % bdfm summary
            elseif model_type == 2
                self.add_bdfm_summary();    
            % midas summary
            elseif model_type == 3
                self.add_midas_summary();                  
            end
        end


        function add_nowcasting_header(self)
            % recover model name and create header
            model_name = self.complementary_information.model_name;
            self.estimation_summary = [self.estimation_summary;cu.model_header(model_name)];
        end


        function add_nowcasting_estimation_header(self)
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
            left_element = sprintf(['Frequency:' sprintf('%+28s', frequency)]);
            right_element = '                                      ';
            lines = [lines; [left_element '    ' right_element]];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_mfbvar_summary(self)
            % make list of regressors
            self.mfbvar_regressors_and_index();
            % loop over equations
            for i=1:self.model.n
                % add coefficient summary
                self.add_mfbvar_coefficient_summary(i);
                % residual and shock variance
                self.add_mfbvar_shock_variance_summary(i);
                % in-sample fit criteria
                self.add_mfbvar_insample_evaluation(i);
            end
            % residual variance-covariance matrix
            self.add_mfbvar_residual_matrix_summary();
            % structural shocks variance-covariance matrix
            self.add_mfbvar_shock_matrix_summary();
            % structural identification matrix
            self.add_mfbvar_structural_identification_matrix_summary();
            % add forecast evaluation criteria, if relevant
            self.add_mfbvar_forecast_evaluation();
        end


        function mfbvar_regressors_and_index(self)
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
            self.mfbvar_regressors = regressors;
            self.mfbvar_coefficient_index = coefficient_index;
        end


        function add_mfbvar_coefficient_summary(self, i)
            lines = string([]);
            endogenous_variables = self.complementary_information.endogenous_variables;
            credibility_level = self.model.credibility_level;
            lines = [lines;cu.equation_header(['Equation: ' char(endogenous_variables(i))])];
            lines = [lines;cu.coefficient_header(credibility_level)];
            lines = [lines;cu.string_line('VAR coefficients beta:')];
            % loop over equation coefficients
            coefficient_index = self.mfbvar_coefficient_index;
            regressors = self.mfbvar_regressors;
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


        function add_mfbvar_shock_variance_summary(self, i)
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


        function add_mfbvar_insample_evaluation(self, i)
            % initiate lines
            lines = string([]);
            % check if in-sample evaluation has been conducted
            if ~isempty(self.model.insample_evaluation)
                lines = [lines;cu.hyphen_dashed_line()];
                ssr = self.model.insample_evaluation.ssr(i);
                r2 = self.model.insample_evaluation.r2(i);
                adj_r2 = self.model.insample_evaluation.adj_r2(i); 
                aic = [];
                bic = [];
                m_y = [];
                lines = [lines;cu.insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic)];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function add_mfbvar_residual_matrix_summary(self)
            Sigma = self.model.Sigma_estimates;
            n = self.model.n;
            endogenous_variables = self.complementary_information.endogenous_variables;  
            lines = string([]);
            lines = [lines;cu.equation_header('Residual variance-covariance Sigma')];
            lines = [lines;cu.variance_covariance_summary(Sigma, n, endogenous_variables, 'var.Sigma_estimates')];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_mfbvar_shock_matrix_summary(self)
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


        function add_mfbvar_structural_identification_matrix_summary(self)
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


        function add_mfbvar_forecast_evaluation(self)
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
                    lines = [lines;cu.forecast_evaluation_line(char(endogenous_variables(i)), ...
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


        function add_bdfm_summary(self)
            % make list of regressors
            self.bdfm_regressors_and_index();
            % add loadings summary
            self.add_bdfm_loadings_summary();
            % add factor summary
            self.add_bdfm_factor_summary();
            % add residual summary
            self.add_bdfm_residual_summary();
            % add shock variance summary
            self.add_bdfm_shock_variance_summary();
            % add forecast evaluation criteria, if relevant
            self.add_bdfm_forecast_evaluation(); 
        end


        function bdfm_regressors_and_index(self)
            k = self.model.k;
            l = self.model.l;
            m = self.model.m;
            q = self.model.q;
            p = self.model.p;
            r = self.model.r;
            [loadings_regressors factor_regressors residual_regressors] = cu.make_dfm_regressors(m, q, p, r);
            [loadings_index factor_index] = cu.make_dfm_index(k, l, m, q, p);
            self.bdfm_loadings_regressors = loadings_regressors;
            self.bdfm_factor_regressors = factor_regressors;
            self.bdfm_residual_regressors = residual_regressors;
            self.bdfm_loadings_index = loadings_index;
            self.bdfm_factor_index = factor_index;
        end


        function add_bdfm_loadings_summary(self)
            lines = string([]);
            endogenous_variables = self.complementary_information.endogenous_variables;
            credibility_level = self.model.credibility_level;
            loadings_regressors = self.bdfm_loadings_regressors;
            loadings_index = self.bdfm_loadings_index;
            if ~isempty(self.model.insample_evaluation)
                evaluation = true;
                insample_evaluation = self.model.insample_evaluation;
            else
                evaluation = false;
            end
            for i=1:self.model.n
                lines = [lines;cu.equation_header(['Loadings equation: ' char(endogenous_variables(i))])];
                lines = [lines;cu.coefficient_header(credibility_level)];
                lines = [lines;cu.string_line('Loadings coefficients lambda:')];
                % loop over equation coefficients
                for j=1:self.model.l
                    regressor = loadings_regressors(j);
                    index = loadings_index(j);
                    coefficient = self.model.lambda_estimates(i,index,1);
                    standard_deviation = self.model.lambda_estimates(i,index,4);
                    lower_bound = self.model.lambda_estimates(i,index,2);
                    upper_bound = self.model.lambda_estimates(i,index,3);
                    lines = [lines;cu.string_line(cu.parameter_estimate_line(regressor,...
                             coefficient, standard_deviation, lower_bound, upper_bound))];
                end
                if evaluation
                    lines = [lines;cu.hyphen_dashed_line()];
                    ssr = self.model.insample_evaluation.ssr(i);
                    r2 = self.model.insample_evaluation.r2(i);
                    adj_r2 = self.model.insample_evaluation.adj_r2(i);   
                    lines = [lines;cu.insample_evaluation_lines(ssr, r2, adj_r2, [], [], [])];
                end
            end
            self.estimation_summary = [self.estimation_summary;lines];
        end

    
        function add_bdfm_factor_summary(self)
            lines = string([]);
            variables = "factor " + string(1:self.model.m);
            credibility_level = self.model.credibility_level;
            factor_regressors = self.bdfm_factor_regressors;
            factor_index = self.bdfm_factor_index;
            for i=1:self.model.m
                lines = [lines;cu.equation_header(['Factor equation: ' char(variables(i))])];
                lines = [lines;cu.coefficient_header(credibility_level)];
                lines = [lines;cu.string_line('VAR coefficients beta:')];  
                % loop over equation coefficients
                for j=1:self.model.k
                    regressor = factor_regressors(j);
                    index = factor_index(j);
                    coefficient = self.model.beta_estimates(index,i,1);
                    standard_deviation = self.model.beta_estimates(index,i,4);
                    lower_bound = self.model.beta_estimates(index,i,2);
                    upper_bound = self.model.beta_estimates(index,i,3);
                    lines = [lines;cu.string_line(cu.parameter_estimate_line(regressor,...
                             coefficient, standard_deviation, lower_bound, upper_bound))];
                end
            end
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_bdfm_residual_summary(self)
            if self.model.r > 0
                lines = string([]);
                endogenous_variables = self.complementary_information.endogenous_variables;
                credibility_level = self.model.credibility_level;
                residual_regressors = self.bdfm_residual_regressors;
                for i=1:self.model.n
                    lines = [lines;cu.equation_header(['Residual equation: ' char(endogenous_variables(i))])];
                    lines = [lines;cu.coefficient_header(credibility_level)];
                    lines = [lines;cu.string_line('AR coefficients gamma:')]; 
                    % loop over equation coefficients
                    for j=1:self.model.r
                        regressor = residual_regressors(j);
                        index = j;
                        coefficient = self.model.gamma_estimates(i,index,1);
                        standard_deviation = self.model.gamma_estimates(i,index,4);
                        lower_bound = self.model.gamma_estimates(i,index,2);
                        upper_bound = self.model.gamma_estimates(i,index,3);
                        lines = [lines;cu.string_line(cu.parameter_estimate_line(regressor,...
                                 coefficient, standard_deviation, lower_bound, upper_bound))];
                    end
                end
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end

    
        function add_bdfm_shock_variance_summary(self)
            lines = string([]);
            lines = [lines;cu.equation_header('Structural shock variance')];
            residual_variance = self.model.sigma;
            shock_variance = self.model.omega;
            lines = [lines;cu.variance_line(residual_variance, shock_variance)];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_bdfm_forecast_evaluation(self)
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
                    lines = [lines;cu.forecast_evaluation_line(char(endogenous_variables(i)), ...
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

    
        function add_midas_summary(self)

            % add constant summary
            self.add_midas_constant_summary();
            % add endogenous lags summary
            self.add_midas_endogenous_summary();  
            % add exogenous lags summary
            self.add_midas_exogenous_summary();          
            % add residual variance summary
            self.add_midas_residual_variance_summary();  
            % add in-sample evaluation criteria
            self.add_midas_insample_evaluation();      
            % add forecast evaluation criteria, if relevant
            self.add_midas_forecast_evaluation();
        end


        function add_midas_constant_summary(self)
            lines = string([]);
            lines = [lines;cu.equation_header('Constant')];
            lines = [lines;cu.coefficient_header(self.model.credibility_level)];
            coefficient = self.model.beta_estimates(1,1);
            standard_deviation = self.model.beta_estimates(1,4);
            lower_bound = self.model.beta_estimates(1,2);
            upper_bound = self.model.beta_estimates(1,3);
            lines = [lines;cu.string_line(cu.parameter_estimate_line('Constant c:',...
                     coefficient, standard_deviation, lower_bound, upper_bound))];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_midas_endogenous_summary(self)
            if self.model.p > 0
                lines = string([]);
                variable = char(self.complementary_information.endogenous_variables(1));
                lines = [lines;cu.equation_header(['Endogenous lags: ' variable])];
                lines = [lines;cu.coefficient_header(self.model.credibility_level)];
                lines = [lines;cu.string_line('Autoregressive coefficients alpha:')];                
                % loop over equation coefficients
                for i=1:self.model.p
                    regressor = [variable ' (-' num2str(i) ')'];
                    coefficient = self.model.beta_estimates(1+i,1);
                    standard_deviation = self.model.beta_estimates(1+i,4);
                    lower_bound = self.model.beta_estimates(1+i,2);
                    upper_bound = self.model.beta_estimates(1+i,3);
                    lines = [lines;cu.string_line(cu.parameter_estimate_line(regressor,...
                             coefficient, standard_deviation, lower_bound, upper_bound))];
                end
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function add_midas_exogenous_summary(self)
            lines = string([]);
            index = self.model.p + 1;
            exogenous_variables = self.complementary_information.exogenous_variables;
            for i=1:self.model.n
                lines = [lines;cu.equation_header(['Exogenous lags: ' char(exogenous_variables(i))])];
                lines = [lines;cu.coefficient_header(self.model.credibility_level)];
                lines = [lines;cu.string_line('Autoregressive coefficients beta:')]; 
                % loop over equation coefficients
                for j=1:self.model.p_(i)+1
                    index = index + 1;
                    if j == 1
                        regressor = char(exogenous_variables(i));
                    else
                        regressor = [char(exogenous_variables(i)) ' (-' num2str(j-1) ')'];
                    end
                    coefficient = self.model.beta_estimates(index,1);
                    standard_deviation = self.model.beta_estimates(index,4);
                    lower_bound = self.model.beta_estimates(index,2);
                    upper_bound = self.model.beta_estimates(index,3);
                    lines = [lines;cu.string_line(cu.parameter_estimate_line(regressor,...
                             coefficient, standard_deviation, lower_bound, upper_bound))];
                end
            end
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_midas_residual_variance_summary(self)
            lines = string([]);
            lines = [lines;cu.equation_header('Residual variance sigma')];
            residual_variance = self.model.sigma_estimates(1);
            shock_variance = ' ';
            lines = [lines;cu.variance_line(residual_variance, shock_variance)];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_midas_insample_evaluation(self)
            lines = string([]);
            lines = [lines;cu.equation_header('In-sample evaluation criteria')];
            if ~isempty(self.model.insample_evaluation)            
                ssr = self.model.insample_evaluation.ssr;
                r2 = self.model.insample_evaluation.r2;
                adj_r2 = self.model.insample_evaluation.adj_r2;
                lines = [lines;[cu.insample_evaluation_lines(ssr, r2, adj_r2, [], [], [])]];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function add_midas_forecast_evaluation(self)
            % initiate lines
            lines = string([]);
            % check if forecast evaluation has been conducted
            if ~isempty(self.model.forecast_evaluation_criteria) 
                rmse = self.model.forecast_evaluation_criteria.rmse;
                mae = self.model.forecast_evaluation_criteria.mae;
                mape = self.model.forecast_evaluation_criteria.mape;
                theil_u = self.model.forecast_evaluation_criteria.theil_u;
                bias = self.model.forecast_evaluation_criteria.bias;
                log_score = self.model.forecast_evaluation_criteria.log_score;
                crps = self.model.forecast_evaluation_criteria.crps;
                lines = [lines;cu.equation_header('Forecast evaluation criteria')];
                lines = [lines;cu.forecast_evaluation_lines(rmse, mae, mape, theil_u, bias, log_score, crps)];
            end
            lines = [lines;cu.equal_dashed_line()];
            self.estimation_summary = [self.estimation_summary;lines];
        end
         

        function make_nowcasting_application_summary(self)
            % get model type, midas regression, mfbvar or bdfm
            model_type = self.complementary_information.model_type;   
            % mfbvar application summary
            if model_type == 1
                self.make_mfbvar_application_summary();
            % bdfm application summary
            elseif model_type == 2
                self.make_bdfm_application_summary();    
            % midas application summary
            elseif model_type == 3
                self.make_midas_application_summary();                  
            end
        end


        function save_nowcasting_application(self, path)
            % get model type, midas regression, mfbvar or bdfm
            model_type = self.complementary_information.model_type; 
            % mfbvar application summary
            if model_type == 1
                self.save_mfbvar_application(path);
            % bdfm summary
            elseif model_type == 2
                self.save_bdfm_application(path);
            % midas summary
            elseif model_type == 3
                self.save_midas_application(path);
            end
        end


        function make_mfbvar_application_summary(self)
            % in-sample fit measures
            self.make_mfbvar_insample_fit_summary();
            % forecasts
            self.make_mfbvar_forecast_summary();        
            % conditional forecasts
            self.make_mfbvar_conditional_forecast_summary();  
            % impulse response function
            self.make_mfbvar_irf_summary();        
            % forecast error variance decomposition
            self.make_mfbvar_fevd_summary();    
            % historical decomposition
            self.make_mfbvar_hd_summary();
        end


        function save_mfbvar_application(self, path)
            % save in-sample fit
            self.save_mfbvar_insample_fit_summary(path);
            % save forecasts
            self.save_mfbvar_forecast_summary(path);   
            % save conditional forecasts
            self.save_mfbvar_conditional_forecast_summary(path);      
            % save impulse response function
            self.save_mfbvar_irf_summary(path);
            % save forecast error variance decomposition
            self.save_mfbvar_fevd_summary(path);
            % save historical decomposition
            self.save_mfbvar_hd_summary(path);
        end


        function make_mfbvar_insample_fit_summary(self)
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


        function make_mfbvar_forecast_summary(self)
            % run only if forecast has been run
            if ~isempty(self.model.forecast_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                p = self.model.p;
                Y = self.model.Y;
                insample_index = self.complementary_information.dates(p+1:end);
                forecasts = self.model.forecast_estimates;
                forecast_index = self.complementary_information.forecast_dates;
                index = [insample_index;forecast_index];
                forecast_table = table(index);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    header = [string([variable '_actual']) [variable '_med'] [variable '_low'] [variable '_upp']];
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


        function make_mfbvar_conditional_forecast_summary(self)
            % run only if conditional forecast has been run
            if isprop(self.model, 'conditional_forecast_estimates') && ~isempty(self.model.conditional_forecast_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                p = self.model.p;
                Y = self.model.Y;
                insample_index = self.complementary_information.dates(p+1:end);
                forecast_index = self.complementary_information.conditional_forecast_dates;
                forecasts = self.model.conditional_forecast_estimates;
                index = [insample_index;forecast_index];
                forecast_table = table(index);
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


        function make_mfbvar_irf_summary(self)
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


        function make_mfbvar_fevd_summary(self)
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


        function make_mfbvar_hd_summary(self)
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


        function save_mfbvar_insample_fit_summary(self, path)
            if isfield(self.application_summary,'insample_fit')
                insample_fit_summary = self.application_summary.insample_fit;
                full_path = fullfile(path, 'insample_fit.csv');
                writetable(insample_fit_summary, full_path);
            end
        end


        function save_mfbvar_forecast_summary(self, path)
            if isfield(self.application_summary,'forecast')
                forecast_summary = self.application_summary.forecast;
                full_path = fullfile(path, 'forecast.csv');
                writetable(forecast_summary, full_path);
            end
        end


        function save_mfbvar_conditional_forecast_summary(self, path)
            if isfield(self.application_summary,'conditional_forecast')
                conditional_forecast_summary = self.application_summary.conditional_forecast;
                full_path = fullfile(path, 'conditional_forecast.csv');
                writetable(conditional_forecast_summary, full_path);
            end
        end


        function save_mfbvar_irf_summary(self, path)
            if isfield(self.application_summary,'irf')
                irf_summary = self.application_summary.irf;
                full_path = fullfile(path, 'irf.csv');
                writetable(irf_summary, full_path);
            end
        end


        function save_mfbvar_fevd_summary(self, path)
            if isfield(self.application_summary,'fevd')
                fevd_summary = self.application_summary.fevd;
                full_path = fullfile(path, 'fevd.csv');
                writetable(fevd_summary, full_path);
            end
        end


        function save_mfbvar_hd_summary(self, path)
            if isfield(self.application_summary,'hd')
                hd_summary = self.application_summary.hd;
                full_path = fullfile(path, 'hd.csv');
                writetable(hd_summary, full_path);
            end
        end


        function make_bdfm_application_summary(self)
            % in-sample fit measures
            self.make_bdfm_insample_fit_summary();
            % forecasts
            self.make_bdfm_forecast_summary();        
            % impulse response function
            self.make_bdfm_irf_summary();        
            % forecast error variance decomposition
            self.make_bdfm_fevd_summary();    
            % historical decomposition
            self.make_bdfm_hd_summary();
        end


        function save_bdfm_application(self, path)
            % save in-sample fit
            self.save_bdfm_insample_fit_summary(path);
            % save forecasts
            self.save_bdfm_forecast_summary(path);   
            % save conditional forecasts
            self.save_bdfm_conditional_forecast_summary(path);      
            % save impulse response function
            self.save_bdfm_irf_summary(path);
            % save forecast error variance decomposition
            self.save_bdfm_fevd_summary(path);
            % save historical decomposition
            self.save_bdfm_hd_summary(path);
        end


        function make_bdfm_insample_fit_summary(self)
            if ~isempty(self.model.fitted_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                fitted = self.model.fitted_estimates;
                residuals = self.model.residual_estimates;
                n = self.model.n;
                index = self.complementary_information.dates;
                fitted_table = table(index);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    header = [string([variable '_fit_med']) [variable '_fit_low'] [variable '_fit_upp'] ...
                              [variable '_res_med'] [variable '_res_low'] [variable '_res_upp']];
                    fit_med = fitted(:,i,1);
                    fit_low = fitted(:,i,2);
                    fit_upp = fitted(:,i,3);
                    res_med = residuals(:,i,1);
                    res_low = residuals(:,i,2);
                    res_upp = residuals(:,i,3); 
                    variable_table = table(fit_med, fit_low, fit_upp, res_med, res_low, res_upp);
                    variable_table.Properties.VariableNames = header;
                    fitted_table = [fitted_table variable_table];
                end
                self.application_summary.insample_fit = fitted_table;
            end
            factors = self.model.f_estimates;
            m = self.model.m;
            factor_table = table(index);
            for i=1:m
                variable = ['factor_'  num2str(i)];
                header = [string([variable '_med']) [variable '_low'] [variable '_upp']];
                factor_med = fitted(:,i,1);
                factor_low = fitted(:,i,2);
                factor_upp = fitted(:,i,3);
                variable_table = table(factor_med, factor_low, factor_upp);
                variable_table.Properties.VariableNames = header;
                factor_table = [factor_table variable_table];
            end
            self.application_summary.factors = factor_table;
        end


        function make_bdfm_forecast_summary(self)
            if ~isempty(self.model.forecast_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                Y = self.model.fitted_estimates(:,:,1);
                insample_index = self.complementary_information.dates;
                forecasts = self.model.forecast_estimates;
                forecast_index = self.complementary_information.forecast_dates;
                index = [insample_index;forecast_index];
                forecast_table = table(index);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    header = [string([variable '_actual']) [variable '_med'] [variable '_low'] [variable '_upp']];
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


        function make_bdfm_irf_summary(self)
            if ~isempty(self.model.irf_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                m = self.model.m;
                irf = permute(self.model.irf_estimates, [3 4 1 2]);
                index = (1:size(irf,1))';
                irf_table = table(index);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    for j=1:m+1
                        if j <= m
                            shock = ['factor' num2str(j) '_shock'];
                        elseif j == m+1
                            shock = 'own_shock';
                        end
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
        end


        function make_bdfm_fevd_summary(self)
            if ~isempty(self.model.fevd_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                m = self.model.m;
                fevd = permute(self.model.fevd_estimates, [3 4 1 2]);
                index = (1:size(fevd,1))';
                fevd_table = table(index);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    for j=1:m+1
                        if j <= m
                            shock = ['factor' num2str(j) '_shock'];
                        elseif j == m+1
                            shock = 'own_shock';
                        end
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


        function make_bdfm_hd_summary(self)
            if ~isempty(self.model.hd_estimates)
                endogenous_variables = self.complementary_information.endogenous_variables;
                n = self.model.n;
                m = self.model.m;
                hd = permute(self.model.hd_estimates, [3 4 1 2]);
                index = self.complementary_information.dates;
                hd_table = table(index);
                for i=1:n
                    variable = char(endogenous_variables(i));
                    for j=1:m+1
                        if j <= m
                            shock = ['factor' num2str(j) '_shock'];
                        elseif j == m+1
                            shock = 'own_shock';
                        end
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


        function save_bdfm_insample_fit_summary(self, path)
            if isfield(self.application_summary,'insample_fit')
                insample_fit_summary = self.application_summary.insample_fit;
                full_path = fullfile(path, 'insample_fit.csv');
                writetable(insample_fit_summary, full_path);
            end
        end


        function save_bdfm_forecast_summary(self, path)
            if isfield(self.application_summary,'forecast')
                forecast_summary = self.application_summary.forecast;
                full_path = fullfile(path, 'forecast.csv');
                writetable(forecast_summary, full_path);
            end
        end


        function save_bdfm_conditional_forecast_summary(self, path)
            if isfield(self.application_summary,'conditional_forecast')
                conditional_forecast_summary = self.application_summary.conditional_forecast;
                full_path = fullfile(path, 'conditional_forecast.csv');
                writetable(conditional_forecast_summary, full_path);
            end
        end


        function save_bdfm_irf_summary(self, path)
            if isfield(self.application_summary,'irf')
                irf_summary = self.application_summary.irf;
                full_path = fullfile(path, 'irf.csv');
                writetable(irf_summary, full_path);
            end
        end


        function save_bdfm_fevd_summary(self, path)
            if isfield(self.application_summary,'fevd')
                fevd_summary = self.application_summary.fevd;
                full_path = fullfile(path, 'fevd.csv');
                writetable(fevd_summary, full_path);
            end
        end


        function save_bdfm_hd_summary(self, path)
            if isfield(self.application_summary,'hd')
                hd_summary = self.application_summary.hd;
                full_path = fullfile(path, 'hd.csv');
                writetable(hd_summary, full_path);
            end
        end


        function make_midas_application_summary(self)
            % in-sample fit measures
            self.make_midas_insample_fit_summary();
            % forecasts
            self.make_midas_forecast_summary();        
        end


        function save_midas_application(self, path)
            % save in-sample fit
            self.save_midas_insample_fit_summary(path);
            % save forecasts
            self.save_midas_forecast_summary(path);   
        end


        function make_midas_insample_fit_summary(self)
            if ~isempty(self.model.fitted_estimates)
                index = self.complementary_information.dates(end-size(self.model.y,1)+1:end);
                fitted_table = table(index);
                actual = self.model.y;
                fitted_med = self.model.fitted_estimates(:,1);
                fitted_low = self.model.fitted_estimates(:,2);
                fitted_upp = self.model.fitted_estimates(:,3);
                residual_med = self.model.residual_estimates(:,1);
                residual_low = self.model.residual_estimates(:,2);
                residual_upp = self.model.residual_estimates(:,3);
                fitted_table = table(index, actual, fitted_med, fitted_low, fitted_upp,...
                                     residual_med, residual_low, residual_upp);
                self.application_summary.insample_fit = fitted_table;
            end
        end


        function make_midas_forecast_summary(self)
            if numel(self.model.forecast_estimates) > 1 || numel(self.model.forecast_estimates{1}) > 0
                forecasts = cell2mat(cellfun(@(x) x.', self.model.forecast_estimates, ...
                            'UniformOutput', false));
                insample_index = self.complementary_information.dates(end-size(self.model.y,1)+1:end);
                forecast_index = self.complementary_information.forecast_dates;
                variable = char(self.complementary_information.endogenous_variables(1));
                header = [string([variable '_actual']) [variable '_med'] [variable '_low'] [variable '_upp']];
                insample_table = table(self.model.y);
                insample_table.Var2(:) = NaN;
                insample_table.Var3(:) = NaN;
                insample_table.Var4(:) = NaN;
                insample_table(end,:) = insample_table(end,1);
                Var2  = forecasts(:,1);
                Var3  = forecasts(:,2);
                Var4  = forecasts(:,3);
                prediction_table = table(Var2, Var3, Var4);
                prediction_table.Var1(:) = NaN;
                variable_table = [insample_table;prediction_table];
                variable_table.Properties.VariableNames = header;
                index = [insample_index;forecast_index];
                forecast_table = table(index);
                forecast_table = [forecast_table variable_table];
                self.application_summary.forecast = forecast_table;
            end
        end


        function save_midas_insample_fit_summary(self, path)
            if isfield(self.application_summary,'insample_fit')
                insample_fit_summary = self.application_summary.insample_fit;
                full_path = fullfile(path, 'insample_fit.csv');
                writetable(insample_fit_summary, full_path);
            end
        end


        function save_midas_forecast_summary(self, path)
            if isfield(self.application_summary,'forecast')
                forecast_summary = self.application_summary.forecast;
                full_path = fullfile(path, 'forecast.csv');
                writetable(forecast_summary, full_path);
            end
        end


    end

end
