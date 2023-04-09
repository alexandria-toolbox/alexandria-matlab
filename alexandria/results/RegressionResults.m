classdef RegressionResults < handle & Results
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        estimation_start
        estimation_complete
        endogenous
        exogenous
        frequency
        start_date
        end_date
        project_path
        data_file
        progress_bar
        create_graphics
        save_results
        regression_type
        iterations
        burnin
        model_credibility
        b
        V
        alpha
        delta
        g
        Q
        tau
        thinning
        thinning_frequency
        Z_variables
        q
        p
        H
        constant
        b_constant
        V_constant
        trend
        b_trend
        V_trend
        quadratic_trend
        b_quadratic_trend
        V_quadratic_trend
        insample_fit
        marginal_likelihood
        hyperparameter_optimization
        optimization_type
        forecast
        forecast_credibility
        forecast_file
        forecast_evaluation
        actual
        insample_dates
        forecast_dates
        k
        n
        beta
        sigma
        h
        gamma
        phi
        fitted
        residuals
        insample_evaluation
        m_y
        estimates_forecasts
        forecast_evaluation_criteria
        summary
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    


        function self = RegressionResults()
            % save estimation start time as attribute
            self.estimation_start = now;
            % display Alexandria header
            self.print_alexandria_header();
            % display Alexandria start message
            self.print_start_message();           
        end
        
        
        function result_summary(self, ip, lr)
            % save estimation end time as attribute
            self.estimation_complete = now;
            % gather information from input processor
            self.input_information(ip);
            % then gather information from regression model
            self.regression_information(lr);
            % print completion message
            self.print_completion_message();
            % build the string list for result summary
            self.generate_result_summary();
            % print and save regression result summary
            self.print_and_save_summary();
        end
        
        
        function settings_summary(self)
            if self.save_results
                % initiate string list
                self.settings = [];
                % add settings header
                self.add_settings_header();
                % add tab 1 settings
                self.add_tab_1_settings();
                % add tab 2 settings
                self.add_tab_2_settings();
                % add tab 3 settings
                self.add_tab_3_settings();
                % print and save regression result summary
                self.save_settings();
            end
        end
        
        
        function application_summary(self)
            if self.save_results
                % actual, fitted and residuals
                if self.insample_fit
                    self.save_fitted_and_residuals();
                end
                % forecasts
                if self.forecast
                    self.save_forecasts();
                end
            end
        end
        

    end
    
    
    methods (Access = protected, Hidden = true)
        
        
        function input_information(self, ip)
            % input processor information: endogenous
            self.endogenous = ip.endogenous_variables;
            % input processor information: exogenous
            self.exogenous = ip.exogenous_variables;
            % input processor information: data frequency
            self.frequency = ip.frequency;
            % input processor information: sample dates
            self.start_date = ip.start_date; 
            self.end_date = ip.end_date;
            % input processor information: project folder
            self.project_path = ip.project_path;
            % input processor information: data file
            self.data_file = ip.data_file;
            % input processor information: progress bar
            self.progress_bar = ip.progress_bar;
            % input processor information: graphics and figures
            self.create_graphics = ip.create_graphics;
            % input processor information: save results
            self.save_results = ip.save_results;
            % input processor information: model
            self.regression_type = ip.regression_type; 
            % input processor information: iterations
            self.iterations = ip.iterations;
            % input processor information: burn-in
            self.burnin = ip.burnin;
            % input processor information: credibility level
            self.model_credibility = ip.model_credibility;
            % input processor information: hyperparameter b
            self.b = ip.b;
            % input processor information: hyperparameter V
            self.V = ip.V;
            % input processor information: hyperparameter alpha
            self.alpha = ip.alpha;
            % input processor information: hyperparameter delta
            self.delta = ip.delta;
            % input processor information: hyperparameter g
            self.g = ip.g;
            % input processor information: hyperparameter Q
            self.Q = ip.Q;   
            % input processor information: hyperparameter tau
            self.tau = ip.tau;       
            % input processor information: thinning
            self.thinning = ip.thinning;
            % input processor information: thinning frequency
            self.thinning_frequency = ip.thinning_frequency;
            % input processor information: Z regressors
            self.Z_variables = ip.Z_variables;
            % input processor information: q
            self.q = ip.q;
            % input processor information: p        
            self.p = ip.p;
            % input processor information: H
            self.H = ip.H;
            % input processor information: constant
            self.constant = ip.constant;
            % input processor information: b constant        
            self.b_constant = ip.b_constant;
            % input processor information: V constant        
            self.V_constant = ip.V_constant;
            % input processor information: trend
            self.trend = ip.trend;
            % input processor information: b trend        
            self.b_trend = ip.b_trend;
            % input processor information: V trend        
            self.V_trend = ip.V_trend;
            % input processor information: quadratic trend
            self.quadratic_trend = ip.quadratic_trend;
            % input processor information: b trend        
            self.b_quadratic_trend = ip.b_quadratic_trend;
            % input processor information: V trend        
            self.V_quadratic_trend = ip.V_quadratic_trend;
            % input processor information: in-sample fit
            self.insample_fit = ip.insample_fit;
            % input processor information: marginal likelihood
            self.marginal_likelihood = ip.marginal_likelihood;
            % input processor information: hyperparameter optimization
            self.hyperparameter_optimization = ip.hyperparameter_optimization;
            % input processor information: optimization type
            self.optimization_type = ip.optimization_type;
            % input processor information: forecasts
            self.forecast = ip.forecast;
            % input processor information: forecasts credibility level
            self.forecast_credibility = ip.forecast_credibility;
            % input processor information: forecast file      
            self.forecast_file = ip.forecast_file;
            % input processor information: forecast evaluation     
            self.forecast_evaluation = ip.forecast_evaluation;
            % input processor information: actual
            self.actual = ip.endogenous;
            % input processor information: in-sample dates
            self.insample_dates = ip.dates;
            % input processor information: forecast dates
            self.forecast_dates = ip.forecast_dates;
        end
        
        
        function regression_information(self, lr)
            % regression information: dimensions
            self.k = lr.k;
            if self.regression_type == 6
                self.n = lr.T;
            else
                self.n = lr.n;
            end
            % regression information: coefficients
            self.beta = lr.estimates_beta;
            if self.regression_type == 1 || self.regression_type == 2
                self.sigma = lr.sigma;
            else
                self.sigma = lr.estimates_sigma;
            end
            if self.regression_type == 5
                self.h = lr.h;
                self.gamma = lr.estimates_gamma;
            end
            if self.regression_type == 6
                self.q = lr.q;
                self.phi = lr.estimates_phi;
            end
            % regression information: in-sample evaluation
            if self.insample_fit
                self.fitted = lr.estimates_fit;
                self.residuals = lr.estimates_residuals;
                self.insample_evaluation = lr.insample_evaluation;
            end
            % regression information: marginal likelihood
            if self.marginal_likelihood && self.regression_type ~= 1
                self.m_y = lr.m_y;
            end
            % regression information: forecasts
            if self.forecast
                self.estimates_forecasts = lr.estimates_forecasts;
            end
            % regression information: forecast evaluation
            if self.forecast && self.forecast_evaluation
                self.forecast_evaluation_criteria = lr.forecast_evaluation_criteria;
            end
        end
        
        
        function generate_result_summary(self)
            % initiate string list
            self.summary = [];
            % add model header, depending on regression type
            self.add_model_header();
            % add estimation header
            self.add_estimation_header();
            % add coefficient header
            self.add_coefficient_header();
            % add beta coefficients
            self.add_beta_coefficients();
            % add sigma coefficient
            self.add_sigma_coefficient();
            % add gamma coefficients, if relevant
            self.add_gamma_coefficient();       
            % add phi coefficients, if relevant
            self.add_phi_coefficient();      
            % add equal line for separation
            self.add_equal_line();
            % add in-sample evaluation criteria, if relevant
            self.add_insample_evaluation();
            % add forecast evaluation criteria, if relevant
            self.add_forecast_evaluation();
        end

        
        function add_model_header(self)
            if self.regression_type == 1
                model = 'Maximum Likelihood Regression';
            elseif self.regression_type == 2
                model = 'Simple Bayesian Regression';          
            elseif self.regression_type == 3
                model = 'Hierarchical Bayesian Regression';         
            elseif self.regression_type == 4
                model = 'Independent Bayesian Regression';
            elseif self.regression_type == 5
                model = 'Heteroscedastic Bayesian Regression';            
            elseif self.regression_type == 6
                model = 'Autocorrelated Bayesian Regression';
            end
            self.summary = [self.summary; cu.model_header(model)];
        end
        

        function add_estimation_header(self)
            estimation_start = self.estimation_start;
            estimation_complete = self.estimation_complete;
            n = self.n;
            endogenous = self.endogenous(1);
            if self.frequency == 1
                frequency = 'cross-sectional/undated';
            elseif self.frequency == 2
                frequency = 'annual';          
            elseif self.frequency == 3
                frequency = 'quarterly'; 
            elseif self.frequency == 4
                frequency = 'monthly';
            elseif self.frequency == 5
                frequency = 'weekly';
            elseif self.frequency == 6
                frequency = 'daily';
            end
            sample = [convertStringsToChars(self.start_date) ' ' ...
                      convertStringsToChars(self.end_date)];
            self.summary = [self.summary; cu.estimation_header(estimation_start, ...
                            estimation_complete, n, endogenous, frequency, sample)];            
            self.summary = [self.summary; cu.equal_dashed_line()];
        end
        

        function add_coefficient_header(self)     
            self.summary = [self.summary; cu.coefficient_header(self.model_credibility)];
        end        
        
        
        function add_beta_coefficients(self)
            lines = [cu.string_line('regression coefficients beta:')];
            regressors = [];
            if self.constant
                regressors = [regressors "constant"];
            end
            if self.trend
                regressors = [regressors "trend"];
            end            
            if self.quadratic_trend
                regressors = [regressors "quadratic_trend"];
            end
            regressors = [regressors self.exogenous];
            for i = 1:self.k
                regressor = regressors(i);
                coefficient = self.beta(i,2);
                standard_deviation = self.beta(i,4);
                lower_bound = self.beta(i,1);
                upper_bound = self.beta(i,3);            
                lines = [lines; cu.parameter_estimate_line(regressor, coefficient, ...
                         standard_deviation, lower_bound, upper_bound)];
            end
            self.summary = [self.summary;lines];
        end          
        
        
        function add_sigma_coefficient(self)
            lines = cu.hyphen_dashed_line();
            lines = [lines; cu.string_line('residual variance sigma:')];
            lines = [lines; ['resid' repmat(' ', 1, 20) ...
                     cu.format_number(self.sigma) repmat(' ', 1, 45)]];
            self.summary = [self.summary;lines];
        end
        
        
        function add_gamma_coefficient(self)
            if self.regression_type == 5
                lines = cu.hyphen_dashed_line();
                lines = [lines; cu.string_line('heteroscedastic coefficients gamma:')];
                for i = 1:self.h
                    lines = [lines; cu.parameter_estimate_line( ...
                            ['Z[' num2str(i) ']'], self.gamma(i,2), ...
                            self.gamma(i,4), self.gamma(i,1),self.gamma(i,3))];
                end
                self.summary = [self.summary;lines];
            end
        end
        
        
        function add_phi_coefficient(self)
            if self.regression_type == 6
                lines = cu.hyphen_dashed_line();
                lines = [lines; cu.string_line('autocorrelation coefficients phi:')];
                for i = 1:self.q
                    lines = [lines; cu.parameter_estimate_line( ...
                            ['resid[-' num2str(i) ']'], self.phi(i,2), ...
                            self.phi(i,4), self.phi(i,1),self.phi(i,3))];
                end
                self.summary = [self.summary;lines];
            end
        end        
        
        
        function add_insample_evaluation(self)
            if self.insample_fit
                ssr = self.insample_evaluation.ssr;
                r2 = self.insample_evaluation.r2;
                adj_r2 = self.insample_evaluation.adj_r2;
                if self.regression_type == 1
                    aic = self.insample_evaluation.aic;
                    bic = self.insample_evaluation.bic;
                    m_y = [];
                else
                    aic = [];
                    bic = [];
                    if self.marginal_likelihood
                        m_y = self.m_y;
                    else
                        m_y = [];
                    end
                end
                lines = cu.insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic);
                self.summary = [self.summary;lines];
                % add equal line for separation
                self.add_equal_line();               
            end
        end 

        
        function add_forecast_evaluation(self)
            if self.forecast && self.forecast_evaluation
                rmse = self.forecast_evaluation_criteria.rmse;
                mae = self.forecast_evaluation_criteria.mae;
                mape = self.forecast_evaluation_criteria.mape;
                theil_u = self.forecast_evaluation_criteria.theil_u;
                bias = self.forecast_evaluation_criteria.bias;
                if self.regression_type ~= 1
                    log_score = self.forecast_evaluation_criteria.log_score;
                    crps = self.forecast_evaluation_criteria.crps;
                else
                    log_score = [];
                    crps = [];
                end
                lines = cu.forecast_evaluation_lines(rmse, mae, mape, theil_u, bias, log_score, crps);
                self.summary = [self.summary;lines];
                % add equal line for separation
                self.add_equal_line();               
            end
        end        
        

        function add_equal_line(self)
            self.summary = [self.summary; cu.equal_dashed_line()];
        end
        
        
        function add_tab_2_settings(self)
            % initiate lines
            lines = [];
            % header for tab 2
            lines = [lines; "Specifications"];
            lines = [lines; "-----------------"];            
            lines = [lines; " "];              
            % recover elements
            regression_type = self.regression_type;
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
            iterations = num2str(self.iterations);
            burnin = num2str(self.burnin);
            model_credibility = num2str(self.model_credibility);
            b = iu.array_to_char(self.b);
            V = iu.array_to_char(self.V);
            alpha = num2str(self.alpha);
            delta = num2str(self.delta);
            g = iu.array_to_char(self.g);
            Q = iu.array_to_char(self.Q);
            tau = num2str(self.tau);
            thinning = cu.bool_to_string(self.thinning);
            thinning_frequency = num2str(self.thinning_frequency);
            Z_variables = iu.array_to_char(self.Z_variables);
            q = num2str(self.q);
            p = iu.array_to_char(self.p);
            H = iu.array_to_char(self.H);
            constant = self.constant;
            constant_string = cu.bool_to_string(constant);
            b_constant = num2str(self.b_constant);
            V_constant = num2str(self.V_constant);
            trend = self.trend;
            trend_string = cu.bool_to_string(trend);      
            b_trend = num2str(self.b_trend);  
            V_trend = num2str(self.V_trend);
            quadratic_trend = self.quadratic_trend;
            quadratic_trend_string = cu.bool_to_string(quadratic_trend);      
            b_quadratic_trend = num2str(self.b_quadratic_trend);       
            V_quadratic_trend = num2str(self.V_quadratic_trend); 
            insample_fit = cu.bool_to_string(self.insample_fit);
            marginal_likelihood = cu.bool_to_string(self.marginal_likelihood);
            hyperparameter_optimization = self.hyperparameter_optimization;
            hyperparameter_optimization_string = cu.bool_to_string(hyperparameter_optimization);
            optimization_type = num2str(self.optimization_type);
            % other elements for tab 2
            lines = [lines; ['regression type: ' model]];  
            if regression_type == 4 || regression_type == 5 || regression_type == 6
                lines = [lines; ['iterations: ' iterations]];  
                lines = [lines; ['burn-in: ' burnin]];
            end
            lines = [lines; ['credibility level: ' model_credibility]];
            if regression_type ~= 1
                lines = [lines; ['b: ' b]];
                lines = [lines; ['V: ' V]];
            end
            if regression_type ~= 1 && regression_type ~= 2
                lines = [lines; ['alpha: ' alpha]];
                lines = [lines; ['delta: ' delta]];
            end
            if regression_type == 5
                lines = [lines; ['g: ' g]];
                lines = [lines; ['Q: ' Q]];
                lines = [lines; ['tau: ' tau]];
                lines = [lines; ['thinning: ' thinning]];
                lines = [lines; ['thinning frequency: ' thinning_frequency]];
                lines = [lines; ['Z variables: ' Z_variables]];
            end
            if regression_type == 6
                lines = [lines; ['q: ' q]];
                lines = [lines; ['p: ' p]];
                lines = [lines; ['H: ' H]];
            end
            lines = [lines; ['constant: ' constant_string]];
            if constant && regression_type ~= 1  
                lines = [lines; ['b (constant): ' b_constant]];
                lines = [lines; ['V (constant): ' V_constant]];
            end
            lines = [lines; ['trend: ' trend_string]];
            if trend && regression_type ~= 1
                lines = [lines; ['b (trend): ' b_trend]];
                lines = [lines; ['V (trend): ' V_trend]]; 
            end
            lines = [lines; ['quadratic trend: ' quadratic_trend_string]];
            if quadratic_trend && regression_type ~= 1
                lines = [lines; ['b (quadratic trend): ' b_quadratic_trend]];
                lines = [lines; ['V (quadratic trend): ' V_quadratic_trend]];
            end
            lines = [lines; ['in-sample fit: ' insample_fit]];        
            if regression_type ~= 1
                lines = [lines; ['marginal likelihood: ' marginal_likelihood]];
            end
            if regression_type == 2 || regression_type == 3
                lines = [lines; ['hyperparameter optimization: ' hyperparameter_optimization_string]];
                if hyperparameter_optimization
                    lines = [lines; ['optimization type: ' optimization_type]];
                end
            end
            lines = [lines; ' ']; 
            lines = [lines; ' '];
            self.settings = [self.settings; lines];
        end
        
        
        function add_tab_3_settings(self)
            % initiate lines
            lines = [];         
            % header for tab 3
            lines = [lines; "Applications"];
            lines = [lines; "-----------------"];            
            lines = [lines; " "];         
            % recover elements
            forecast = self.forecast;
            forecast_string = cu.bool_to_string(forecast);
            forecast_credibility = num2str(self.forecast_credibility);
            forecast_file = self.forecast_file;
            forecast_evaluation = cu.bool_to_string(self.forecast_evaluation);
            % other elements for tab 3
            lines = [lines; ['forecast: ' forecast_string]];
            if forecast
                lines = [lines; ['forecast credibility: ' forecast_credibility]];
                lines = [lines; ['forecast file: ' forecast_file]];
                lines = [lines; ['forecast evaluation: ' forecast_evaluation]];
            end
            % other elements for tab 3
            self.settings = [self.settings; lines];
        end
        
        
        function save_fitted_and_residuals(self)
            % create index, column labels and data
            index = self.insample_dates;
            actual = self.actual;
            fitted = self.fitted;
            residuals = self.residuals;
            fitted_table = table(index, actual, fitted, residuals);
            % create path to file and save
            fitted_file_path = fullfile(self.project_path, 'results', 'fitted_and_residuals.csv');
            writetable(fitted_table, fitted_file_path);
        end
        
        
        function save_forecasts(self)
            % create index, column labels and data
            index = self.forecast_dates;
            lower_bound = self.estimates_forecasts(:,1);
            median = self.estimates_forecasts(:,2);
            upper_bound = self.estimates_forecasts(:,3);
            forecast_table = table(index, lower_bound, median, upper_bound);
            % create path to file and save
            forecast_file_path = fullfile(self.project_path, 'results', 'forecasts.csv');
            writetable(forecast_table, forecast_file_path);
        end
          

    end
    

end
