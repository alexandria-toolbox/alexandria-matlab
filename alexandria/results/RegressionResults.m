classdef RegressionResults < handle
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    


        function self = RegressionResults()           
        end


    end
    
    
    methods (Access = protected, Hidden = true)


        function complete_regression_information(self)
            % endogenous and exogenous variables
            if ~isfield(self.complementary_information, 'endogenous_variables')
                self.complementary_information.endogenous_variables = ["y"];
            end
            if ~isfield(self.complementary_information, 'exogenous_variables')
                n_exo = size(self.model.exogenous,2);
                self.complementary_information.exogenous_variables = "x"+(1:n_exo);
            end
            % sample dates
            if ~isfield(self.complementary_information, 'dates')
                n = self.model.n;
                self.complementary_information.dates = string(1:n);
            end
            % heteroscedastic variables
            model_type = self.complementary_information.model_type;
            if model_type == 5 && ~isfield(self.complementary_information, 'heteroscedastic_variables')
                self.complementary_information.heteroscedastic_variables = "z"+(1:size(self.model.Z,2));
            end
            % regression options
            if ~isfield(self.complementary_information, 'insample_fit')
                self.complementary_information.insample_fit = '_';
            end
            if ~isfield(self.complementary_information, 'marginal_likelihood')
                self.complementary_information.marginal_likelihood = '_';
            end
            if ~isfield(self.complementary_information, 'hyperparameter_optimization')
                self.complementary_information.hyperparameter_optimization = '_';
            end
            if ~isfield(self.complementary_information, 'optimization_type')
                self.complementary_information.optimization_type = '_';
            end
        end


        function add_regression_tab_2_inputs(self)
            % initiate lines
            lines = string([]);
            % header for tab 2
            lines = [lines;'Specification'];
            lines = [lines;'-----------------'];
            lines = [lines;' '];
            % regression type
            regression_type = self.complementary_information.model_type;
            if regression_type == 1
                model = 'Maximum Likelihood Regression';
            elseif regression_type == 2
                model = 'Simple Bayesian Regression';      
            elseif regression_type == 3
                model = 'Hierarchical Bayesian Regression';         
            elseif regression_type == 4
                model = 'Independent Bayesian Regression';
            elseif regression_type == 5
                model = 'Heteroscedastic Bayesian Regression';            
            elseif regression_type == 6
                model = 'Autocorrelated Bayesian Regression';
            end
            lines = [lines;'regression type: ' model];
            % burn-in and iterations
            if regression_type == 4 || regression_type == 5 || regression_type == 6
                iterations = num2str(self.model.iterations);
                burn = num2str(self.model.burn);
                lines = [lines;'iterations: ' iterations];
                lines = [lines;'burn-in: ' burn];
            end
            % credibility level
            model_credibility = num2str(self.model.credibility_level);
            lines = [lines;'credibility level: ' model_credibility];
            % hyperparameters: b and V
            if regression_type ~= 1
                b = iu.array_to_char(self.model.b);
                V = iu.array_to_char(diag(self.model.V));
                lines = [lines;'b: ' b];
                lines = [lines;'V: ' V];
            end
            % hyperparameters: alpha and delta
            if regression_type ~= 1 && regression_type ~= 2
                alpha = num2str(self.model.alpha);
                delta = num2str(self.model.delta);
                lines = [lines;'alpha: ' alpha];
                lines = [lines;'delta: ' delta];
            end
            % hyperparameters: heteroscedastic regression
            if regression_type == 5      
                g = iu.array_to_char(self.model.g);
                Q = iu.array_to_char(diag(self.model.Q));
                tau = num2str(self.model.tau);
                thinning = cu.bool_to_string(self.model.thinning);
                thinning_frequency = num2str(self.model.thinning_frequency);
                Z_variables = iu.array_to_char(self.complementary_information.heteroscedastic_variables);
                lines = [lines;'g: ' g];
                lines = [lines;'Q: ' Q];    
                lines = [lines;'tau: ' tau];
                lines = [lines;'thinning: ' thinning];
                lines = [lines;'thinning frequency: ' thinning_frequency];
                lines = [lines;'Z variables: ' Z_variables];
            end
            % hyperparameters: autocorrelated regression
            if regression_type == 6     
                q = num2str(self.model.q);
                p = iu.array_to_char(self.model.p);
                H = iu.array_to_char(diag(self.model.H));
                lines = [lines;'q: ' q];
                lines = [lines;'p: ' p];       
                lines = [lines;'H: ' H]; 
            end
            % constant and constant prior  
            constant = cu.bool_to_string(self.model.constant);         
            lines = [lines;'constant: ' constant];
            if regression_type ~= 1
                b_constant = num2str(self.model.b_constant);
                V_constant = num2str(self.model.V_constant); 
                lines = [lines;'b (constant): ' b_constant];           
                lines = [lines;'V (constant): ' V_constant];
            end
            % trend and trend prior  
            trend = cu.bool_to_string(self.model.trend);              
            lines = [lines;'trend: ' trend];
            if regression_type ~= 1
                b_trend = num2str(self.model.b_trend);
                V_trend = num2str(self.model.V_trend); 
                lines = [lines;'b (trend): ' b_trend];          
                lines = [lines;'V (trend): ' V_trend]; 
            end
            % quadratic trend and quadratic trend prior  
            quadratic_trend = cu.bool_to_string(self.model.quadratic_trend);            
            lines = [lines;'quadratic trend: ' quadratic_trend];
            if regression_type ~= 1
                b_quadratic_trend = num2str(self.model.b_quadratic_trend);
                V_quadratic_trend = num2str(self.model.V_quadratic_trend); 
                lines = [lines;'b (quadratic trend): ' b_quadratic_trend];           
                lines = [lines;'V (quadratic trend): ' V_quadratic_trend];
            end
            % in-sample fit
            if islogical(self.complementary_information.insample_fit)
                insample_fit = cu.bool_to_string(self.complementary_information.insample_fit);   
            else
                insample_fit = self.complementary_information.insample_fit;
            end
            lines = [lines;'in-sample fit: ' insample_fit];           
            % marginal likelihood
            if regression_type ~= 1
                if islogical(self.complementary_information.marginal_likelihood)
                    marginal_likelihood = cu.bool_to_string(self.complementary_information.marginal_likelihood);   
                else
                    marginal_likelihood = self.complementary_information.marginal_likelihood;
                end
                lines = [lines;'marginal likelihood: ' marginal_likelihood];
            end
            % hyperparameter optimization
            if regression_type == 2 || regression_type == 3 
                if islogical(self.complementary_information.hyperparameter_optimization)
                    hyperparameter_optimization = cu.bool_to_string(self.complementary_information.hyperparameter_optimization);   
                else
                    hyperparameter_optimization = self.complementary_information.hyperparameter_optimization;
                end
                lines = [lines;'hyperparameter optimization: ' hyperparameter_optimization];
                % optimization type
                optimization_type = self.complementary_information.optimization_type;
                lines = [lines;'optimization type: ' optimization_type];
            end
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end


        function add_regression_tab_3_inputs(self)
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
            % forecast options
            forecast_file = self.complementary_information.forecast_file;
            lines = [lines;'forecast file: ' forecast_file];
            if islogical(self.complementary_information.forecast_evaluation)
                forecast_evaluation = cu.bool_to_string(self.complementary_information.forecast_evaluation);
            else
                forecast_evaluation = self.complementary_information.forecast_evaluation;
            end
            lines = [lines;'forecast evaluation: ' forecast_evaluation];
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end

  
        function make_regression_summary(self)
            % initiate string list
            self.estimation_summary = string([]);
            % add model header
            self.add_regression_header();    
            % add estimation header
            self.add_regression_estimation_header();
            % add coefficient header
            self.add_regression_coefficient_summary();
            % add heteroscedastic coefficients, if relevant
            self.add_heteroscedastic_coefficient(); 
            % add autocorrelated coefficients, if relevant
            self.add_autocorrelated_coefficient();         
            % add sigma coefficient
            self.add_sigma_coefficient();
            % add in-sample evaluation criteria, if relevant
            self.add_regression_insample_evaluation(); 
            % add forecast evaluation criteria, if relevant
            self.add_regression_forecast_evaluation();
        end


        function make_regression_application_summary(self)
            % in-sample fit measures
            self.make_regression_insample_fit_summary();
            % forecasts
            self.make_regression_forecast_summary();
        end


        function save_regression_application(self, path)
            % save in-sample fit
            self.save_regression_insample_fit_summary(path);
            % save forecasts
            self.save_regression_forecast_summary(path);
        end


        function add_regression_header(self)
            % recover model name and create header
            model_name = self.complementary_information.model_name;
            self.estimation_summary = [self.estimation_summary;cu.model_header(model_name)];
        end


        function add_regression_estimation_header(self)
            % initiate lines
            lines = string([]);
            % first row: dependent variable and frequency
            endogenous_variables = convertStringsToChars(self.complementary_information.endogenous_variables(1));
            frequency = self.complementary_information.frequency;
            left_element = sprintf(['Dep. variable:' sprintf('%+24s', cu.shorten_string(endogenous_variables, 20))]);
            right_element = sprintf(['Frequency:' sprintf('%+28s', frequency)]);
            lines = [lines; [left_element '    ' right_element]];
            % second row: estimation sample and estimation start
            sample_start = self.complementary_information.sample_start;
            sample_end = self.complementary_information.sample_end;
            if numel(sample_start) == 0 || numel(sample_end) == 0
                sample = '—';
            else
                sample = [sample_start '  ' sample_end];
            end
            estimation_start = self.complementary_information.estimation_start;
            left_element = sprintf(['Sample:' sprintf('%+31s',sample)]);
            right_element = sprintf(['Est. start:' sprintf('%+27s', estimation_start)]);
            lines = [lines; [left_element '    ' right_element]];    
            % third row: observations and estimation complete   
            n = num2str(self.model.n);
            estimation_end = self.complementary_information.estimation_end;
            left_element = sprintf(['No. observations:' sprintf('%+21s', n)]);
            right_element = sprintf(['Est. complete:' sprintf('%+24s', estimation_end)]);
            lines = [lines; [left_element '    ' right_element]];    
            % equal dashed line
            lines = [lines;cu.equal_dashed_line()];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_regression_coefficient_summary(self)
            % initiate lines
            lines = string([]);
            % coefficients header
            credibility_level = self.model.credibility_level;
            lines = [lines;cu.coefficient_header(credibility_level)];
            lines = [lines;cu.string_line('regression coefficients beta:')];
            % coefficient summary, coefficient by coefficient
            regressors = string([]);
            if self.model.constant
                regressors = [regressors 'constant'];
            end
            if self.model.trend
                regressors = [regressors 'trend'];
            end
            if self.model.quadratic_trend
                regressors = [regressors 'quadratic trend'];
            end
            regressors = [regressors self.complementary_information.exogenous_variables];
            for i=1:self.model.k
                regressor = regressors(i);
                coefficient = self.model.beta_estimates(i,1);
                standard_deviation = self.model.beta_estimates(i,4);
                lower_bound = self.model.beta_estimates(i,2);
                upper_bound = self.model.beta_estimates(i,3);
                lines = [lines; cu.parameter_estimate_line(regressor, ...
                        coefficient, standard_deviation, lower_bound, upper_bound)];
            end
            lines = [lines;cu.hyphen_dashed_line()];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_heteroscedastic_coefficient(self)    
            model_type = self.complementary_information.model_type;
            if model_type == 5
                lines = string([]);
                heteroscedastic_variables = self.complementary_information.heteroscedastic_variables;
                gamma_estimates = self.model.gamma_estimates;
                lines = [lines;cu.string_line('heteroscedastic coefficients gamma:')];
                for i=1:self.model.h
                    lines = [lines;cu.parameter_estimate_line(heteroscedastic_variables(i), ...
                    gamma_estimates(i,1), gamma_estimates(i,4), gamma_estimates(i,2), gamma_estimates(i,3))];
                end
                lines = [lines;cu.hyphen_dashed_line()];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function add_autocorrelated_coefficient(self)
            model_type = self.complementary_information.model_type;
            if model_type == 6
                lines = string([]);
                phi_estimates = self.model.phi_estimates;
                lines = [lines;cu.string_line('autocorrelation coefficients phi:')];          
                for i=1:self.model.q
                    lines = [lines;cu.parameter_estimate_line(['resid[-' num2str(i) ']'], ...
                    phi_estimates(i,1), phi_estimates(i,4), phi_estimates(i,2), phi_estimates(i,3))];
                end
                lines = [lines;cu.hyphen_dashed_line()];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end
   

        function add_sigma_coefficient(self)
            % initiate lines
            lines = string([]);
            % coefficients header
            lines = [lines;cu.string_line('residual variance sigma:')];
            % coefficient summary
            model_type = self.complementary_information.model_type;
            if model_type == 1 || model_type == 2
                sigma = self.model.sigma;
            else
                sigma = self.model.sigma_estimates;
            end
            lines = [lines;['resid' repmat(' ', 1, 20) cu.format_number(sigma) repmat(' ', 1, 45)]];
            lines = [lines;cu.equal_dashed_line()];
            self.estimation_summary = [self.estimation_summary;lines];
        end


        function add_regression_insample_evaluation(self)
            % initiate lines
            lines = string([]);
            % check if in-sample evaluation has been conducted
            if ~isempty(self.model.insample_evaluation)
                ssr = self.model.insample_evaluation.ssr;
                r2 = self.model.insample_evaluation.r2;
                adj_r2 = self.model.insample_evaluation.adj_r2;
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
                lines = [lines;[cu.insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic)]];
                lines = [lines;cu.equal_dashed_line()];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function add_regression_forecast_evaluation(self)
            % initiate lines
            lines = string([]);
            % check if forecast evaluation has been conducted
            if ~isempty(self.model.forecast_evaluation_criteria)
                rmse = self.model.forecast_evaluation_criteria.rmse;
                mae = self.model.forecast_evaluation_criteria.mae;
                mape = self.model.forecast_evaluation_criteria.mape;
                theil_u = self.model.forecast_evaluation_criteria.theil_u;
                bias = self.model.forecast_evaluation_criteria.bias;
                model_type = self.complementary_information.model_type;
                if model_type ~= 1
                    log_score = self.model.forecast_evaluation_criteria.log_score;
                    crps = self.model.forecast_evaluation_criteria.crps;
                else
                    log_score = [];
                    crps = [];
                end
                lines = [lines;cu.forecast_evaluation_lines(rmse, mae, mape, theil_u, bias, log_score, crps)];
                lines = [lines;cu.equal_dashed_line()];
                self.estimation_summary = [self.estimation_summary;lines];
            end
        end


        function make_regression_insample_fit_summary(self)
            % run only if in-sample fit has been run
            if ~isempty(self.model.fitted_estimates)
                % create index, column labels and data
                index = self.complementary_information.dates;
                actual = self.model.endogenous;
                fitted = self.model.fitted_estimates;
                residuals = self.model.residual_estimates;
                fitted_table = table(index, actual, fitted, residuals);
                self.application_summary.insample_fit = fitted_table;
            end
        end


        function make_regression_forecast_summary(self)
            % run only if forecast has been run
            if ~isempty(self.model.forecast_estimates)       
            % create index, column labels and data
                index = (1:size(self.model.forecast_estimates,1))';
                median = self.model.forecast_estimates(:,1);
                lower_bound = self.model.forecast_estimates(:,2);
                upper_bound = self.model.forecast_estimates(:,3);
                forecast_table = table(index, median, lower_bound, upper_bound);
                self.application_summary.forecast = forecast_table;
            end
        end


        function save_regression_insample_fit_summary(self, path)
            if isfield(self.application_summary,'insample_fit')
                insample_fit_summary = self.application_summary.insample_fit;
                full_path = fullfile(path, 'insample_fit.csv'); 
                writetable(insample_fit_summary, full_path);
            end
        end
                
    
        function save_regression_forecast_summary(self, path)
            if isfield(self.application_summary,'forecast')
                forecast_summary = self.application_summary.forecast;
                full_path = fullfile(path, 'forecast.csv');
                writetable(forecast_summary, full_path);
            end
        end

    end

end
