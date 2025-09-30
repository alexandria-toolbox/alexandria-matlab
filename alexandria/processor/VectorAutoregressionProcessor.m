classdef VectorAutoregressionProcessor < handle

    
    properties (GetAccess = public, SetAccess= protected)
        var_type
        var_iterations
        var_burnin
        var_model_credibility
        var_constant
        var_trend
        var_quadratic_trend
        lags
        ar_coefficients
        pi1
        pi2
        pi3
        pi4
        pi5
        pi6
        pi7
        proxy_variables
        lamda
        proxy_prior
        var_insample_fit
        constrained_coefficients
        sums_of_coefficients
        initial_observation
        long_run
        stationary
        var_marginal_likelihood
        var_hyperparameter_optimization
        coefficients_file
        long_run_file
        var_endogenous
        var_exogenous
        var_dates
        proxys
        constrained_coefficients_table
        long_run_table
        var_Z_p
        var_Y_p
        var_forecast_dates
        condition_table
        shock_table
        restriction_table
    end     
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------     
    
    
    methods (Access = public)
        
        
        function self = VectorAutoregressionProcessor()
        end        
        
    
        function vector_autoregression_inputs(self)
            % recover vector autoregression type
            self.var_type = self.get_vector_autoregression_type();
            % recover iterations
            self.var_iterations = self.get_var_iterations();
            % recover burn-in
            self.var_burnin = self.get_var_burnin();
            % recover credibility level for model estimates
            self.var_model_credibility = self.get_var_model_credibility();        
            % recover constant
            self.var_constant = self.get_var_constant(); 
            % recover trend
            self.var_trend = self.get_var_trend(); 
            % recover quadratic trend
            self.var_quadratic_trend = self.get_var_quadratic_trend();         
            % recover lags
            self.lags = self.get_lags();  
            % recover AR coefficients
            self.ar_coefficients = self.get_ar_coefficients();          
            % recover pi1
            self.pi1 = self.get_pi1();  
            % recover pi2
            self.pi2 = self.get_pi2(); 
            % recover pi3
            self.pi3 = self.get_pi3(); 
            % recover pi4
            self.pi4 = self.get_pi4(); 
            % recover pi5
            self.pi5 = self.get_pi5(); 
            % recover pi6
            self.pi6 = self.get_pi6(); 
            % recover pi7
            self.pi7 = self.get_pi7(); 
            % recover proxys
            self.proxy_variables = self.get_proxy_variables();
            % recover lamda
            self.lamda = self.get_lamda(); 
            % recover proxy prior type
            self.proxy_prior = self.get_proxy_prior();
            % recover in-sample fit
            self.var_insample_fit = self.get_var_insample_fit();            
            % recover constrained_coefficients
            self.constrained_coefficients = self.get_constrained_coefficients();
            % recover sums-of-coefficients
            self.sums_of_coefficients = self.get_sums_of_coefficients();     
            % recover dummy initial observation
            self.initial_observation = self.get_initial_observation();
            % recover long-run prior
            self.long_run = self.get_long_run();
            % recover stationary prior
            self.stationary = self.get_stationary();
            % recover marginal likelihood
            self.var_marginal_likelihood = self.get_var_marginal_likelihood();
            % recover hyperparameter optimization
            self.var_hyperparameter_optimization = self.get_var_hyperparameter_optimization();
            % recover constrained coefficients file
            self.coefficients_file = self.get_coefficients_file();
            % recover long-run prior file
            self.long_run_file = self.get_long_run_file();
        end
        

        function vector_autoregression_data(self)
            % print loading message
            if self.progress_bar
                cu.print_message_to_complete('Data loading:');
            end
            % recover in-sample endogenous and exogenous
            [self.var_endogenous, self.var_exogenous, self.var_dates] = self.get_var_insample_data();
            % recover proxy variables
            self.proxys = self.get_proxy_data();
            % recover constrained coefficients
            self.constrained_coefficients_table = self.get_constrained_coefficients_table();
            % recover long run prior
            self.long_run_table = self.get_long_run_table();
            % recover forecast data
            [self.var_Z_p, self.var_Y_p, self.var_forecast_dates] = self.get_var_forecast_data();
            % recover conditional forecast data
            [self.condition_table, self.shock_table] = self.get_condition_table();
            % recover sign restrictions data
            [self.restriction_table] = self.get_restriction_table();
            if self.progress_bar
                cu.print_message('  —  done');
            end
        end


        function make_var_information(self)
            % get sample dates
            self.results_information.dates = self.var_dates;
            % get forecast dates
            self.results_information.forecast_dates = self.var_forecast_dates;
            self.results_information.conditional_forecast_dates = self.var_forecast_dates;
            % get proxy variables for proxy SVAR
            if self.var_type == 7
                self.results_information.proxy_variables = self.proxy_variables;
            else
                self.results_information.proxy_variables = [];
            end
            % get var option: in-sample fit
            self.results_information.insample_fit = self.var_insample_fit;
            % get var option: constrained coefficients
            self.results_information.constrained_coefficients = self.constrained_coefficients;      
            % get var option: sums-of-coefficients
            self.results_information.sums_of_coefficients = self.sums_of_coefficients;
            % get var option: dummy initial observation
            self.results_information.initial_observation = self.initial_observation;
            % get var option: long-run prior
            self.results_information.long_run = self.long_run;
            % get var option: stationary prior
            self.results_information.stationary = self.stationary;
            % get var option: marginal likelihood
            self.results_information.marginal_likelihood = self.var_marginal_likelihood;
            % get var option: hyperparameter optimization
            self.results_information.hyperparameter_optimization = self.var_hyperparameter_optimization;
            % get constrained coefficients file
            self.results_information.coefficients_file = self.coefficients_file;
            % get long run prior file
            self.results_information.long_run_file = self.long_run_file;
        end


        function make_var_graphics_information(self)
            % get sample dates
            self.graphics_information.dates = self.var_dates;
            % get forecast dates
            self.graphics_information.forecast_dates = self.var_forecast_dates;
            self.graphics_information.conditional_forecast_dates = self.var_forecast_dates;
            % get actual data for forecast evaluation, if available
            self.graphics_information.Y_p = self.var_Y_p;
        end


    end
        

    %---------------------------------------------------
    % Methods (Access = private)
    %--------------------------------------------------- 
        
  
    methods (Access = protected, Hidden = true) 
        
        
        function [var_type] = get_vector_autoregression_type(self)
            var_type = self.user_inputs.tab_2_var.var_type;
            if ~ismember(var_type, [1 2 3 4 5 6 7])
                error(['Value error for vector autoregression type. Should be integer between 1 and 7.']);
            end
        end
        

        function [iterations] = get_var_iterations(self)
            iterations = self.user_inputs.tab_2_var.iterations;
            if ~(ischar(iterations) || iu.is_integer(iterations))
                error(['Type error for iterations. Should be integer.']);
            end
            if ~iu.is_empty(iterations) && ischar(iterations)
                if iu.is_digit(iterations)
                    iterations = str2double(iterations);
                else
                    error(['Type error for iterations. Should be positive integer.']);
                end
            end  
            if iu.is_integer(iterations) && iterations <= 0
                error(['Value error for iterations. Should be positive integer.']);
            end            
        end        


        function [burnin] = get_var_burnin(self)
            burnin = self.user_inputs.tab_2_var.burnin;
            if ~(ischar(burnin) || iu.is_integer(burnin))
                error(['Type error for burn-in. Should be integer.']);
            end
            if ~iu.is_empty(burnin) && ischar(burnin)
                if iu.is_digit(burnin)
                    burnin = str2double(burnin);
                else
                    error(['Type error for burn-in. Should be positive integer.']);
                end
            end  
            if iu.is_integer(burnin) && burnin <= 0
                error(['Value error for burn-in. Should be positive integer.']);
            end            
        end 
        

        function [model_credibility] =  get_var_model_credibility(self)
            model_credibility = self.user_inputs.tab_2_var.model_credibility;
            if ~ismember(class(model_credibility), ["char" "double"])
                error(['Type error for model credibility level. Should be float between 0 and 1.']);
            end
            if ischar(model_credibility)
                if isnan(str2double(model_credibility))
                    error(['Type error for model credibility level. Should be float between 0 and 1.']);
                else
                    model_credibility = str2double(model_credibility);
                end
            end
            if model_credibility <= 0 || model_credibility >= 1
                error(['Value error for model credibility level. Should be float between 0 and 1 (not included).']);
            end
        end  
        

        function [constant] = get_var_constant(self)
            constant = self.user_inputs.tab_2_var.constant;
            if ~islogical(constant)
                error(['Type error for constant. Should be boolean.']);
            end
        end         
        

        function [trend] = get_var_trend(self)
            trend = self.user_inputs.tab_2_var.trend;
            if ~islogical(trend)
                error(['Type error for trend. Should be boolean.']);
            end
        end         
        

        function [quadratic_trend] = get_var_quadratic_trend(self)
            quadratic_trend = self.user_inputs.tab_2_var.quadratic_trend;
            if ~islogical(quadratic_trend)
                error(['Type error for quadratic trend. Should be boolean.']);
            end
        end         
        

        function [lags] = get_lags(self)
            lags = self.user_inputs.tab_2_var.lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags <= 0
                error(['Value error for lags. Should be positive integer.']);
            end            
        end          
        

        function [ar_coefficients] = get_ar_coefficients(self)
            ar_coefficients = self.user_inputs.tab_2_var.ar_coefficients;
            if ~ismember(class(ar_coefficients), ["char" "double"])
                error(['Type error for AR coefficients. Should be float or scalar array.']);
            end
            if ischar(ar_coefficients)
                ar_coefficients = iu.char_to_array(ar_coefficients);
                if any(isnan(str2double(ar_coefficients)))
                    error(['Type error for AR coefficients. All elements should be scalars.']);
                else
                    ar_coefficients = str2double(ar_coefficients)';
                end
            end
            if numel(ar_coefficients) ~= numel(self.endogenous_variables) && numel(ar_coefficients) ~= 1
                error(['Dimension error for AR coefficients. Dimension of AR coefficients and endogenous don''t match.']);
            end
            if any(isnan(ar_coefficients)) || any(isinf(ar_coefficients))
                error(['Type error for AR coefficients. All elements should be scalars.']);
            end
        end
        

        function [pi1] = get_pi1(self)
            pi1 = self.user_inputs.tab_2_var.pi1;
            if ~ismember(class(pi1), ["char" "double"])
                error(['Type error for pi1. Should be float or integer.']);
            end
            if ischar(pi1)
                if isnan(str2double(pi1))
                    error(['Type error for pi1. Should be float or integer.']);
                else
                    pi1 = str2double(pi1);
                end
            end
            if pi1 <= 0
                error(['Value error for pi1. Should be strictly positive.']);
            end
        end   
        

        function [pi2] = get_pi2(self)
            pi2 = self.user_inputs.tab_2_var.pi2;
            if ~ismember(class(pi2), ["char" "double"])
                error(['Type error for pi2. Should be float or integer.']);
            end
            if ischar(pi2)
                if isnan(str2double(pi2))
                    error(['Type error for pi2. Should be float or integer.']);
                else
                    pi2 = str2double(pi2);
                end
            end
            if pi2 <= 0
                error(['Value error for pi2. Should be strictly positive.']);
            end
        end  
        

        function [pi3] = get_pi3(self)
            pi3 = self.user_inputs.tab_2_var.pi3;
            if ~ismember(class(pi3), ["char" "double"])
                error(['Type error for pi3. Should be float or integer.']);
            end
            if ischar(pi3)
                if isnan(str2double(pi3))
                    error(['Type error for pi3. Should be float or integer.']);
                else
                    pi3 = str2double(pi3);
                end
            end
            if pi3 <= 0
                error(['Value error for pi3. Should be strictly positive.']);
            end
        end          
        

        function [pi4] = get_pi4(self)
            pi4 = self.user_inputs.tab_2_var.pi4;
            if ~ismember(class(pi4), ["char" "double"])
                error(['Type error for pi4. Should be float or integer.']);
            end
            if ischar(pi4)
                if isnan(str2double(pi4))
                    error(['Type error for pi4. Should be float or integer.']);
                else
                    pi4 = str2double(pi4);
                end
            end
            if pi4 <= 0
                error(['Value error for pi4. Should be strictly positive.']);
            end
        end        
        

        function [pi5] = get_pi5(self)
            pi5 = self.user_inputs.tab_2_var.pi5;
            if ~ismember(class(pi5), ["char" "double"])
                error(['Type error for pi5. Should be float or integer.']);
            end
            if ischar(pi5)
                if isnan(str2double(pi5))
                    error(['Type error for pi5. Should be float or integer.']);
                else
                    pi5 = str2double(pi5);
                end
            end
            if pi5 <= 0
                error(['Value error for pi5. Should be strictly positive.']);
            end
        end        
        

        function [pi6] = get_pi6(self)
            pi6 = self.user_inputs.tab_2_var.pi6;
            if ~ismember(class(pi6), ["char" "double"])
                error(['Type error for pi6. Should be float or integer.']);
            end
            if ischar(pi6)
                if isnan(str2double(pi6))
                    error(['Type error for pi6. Should be float or integer.']);
                else
                    pi6 = str2double(pi6);
                end
            end
            if pi6 <= 0
                error(['Value error for pi6. Should be strictly positive.']);
            end
        end
        

        function [pi7] = get_pi7(self)
            pi7 = self.user_inputs.tab_2_var.pi7;
            if ~ismember(class(pi7), ["char" "double"])
                error(['Type error for pi7. Should be float or integer.']);
            end
            if ischar(pi7)
                if isnan(str2double(pi7))
                    error(['Type error for pi7. Should be float or integer.']);
                else
                    pi7 = str2double(pi7);
                end
            end
            if pi7 <= 0
                error(['Value error for pi7. Should be strictly positive.']);
            end
        end        
        

        function [proxy_variables] = get_proxy_variables(self)
            proxy_variables = self.user_inputs.tab_2_var.proxy_variables;
            if self.var_type == 7
                if iu.is_empty(proxy_variables) || ...
                   ~ismember(class(proxy_variables), ["string" "char"])
                    error(['Type error for proxys. Should non-empty char or string.']);
                end
                proxy_variables = iu.char_to_array(proxy_variables);
            end
        end        
        

        function [lamda] = get_lamda(self)
            lamda = self.user_inputs.tab_2_var.lamda;
            if ~ismember(class(lamda), ["char" "double"])
                error(['Type error for lambda. Should be float or integer.']);
            end
            if ischar(lamda)
                if isnan(str2double(lamda))
                    error(['Type error for lambda. Should be float or integer.']);
                else
                    lamda = str2double(lamda);
                end
            end
            if lamda <= 0 || lamda > 1
                error(['Value error for lambda. Should be float between 0 and 1..']);
            end
        end         
        

        function [proxy_prior] = get_proxy_prior(self)
            proxy_prior = self.user_inputs.tab_2_var.proxy_prior;
            if ~ismember(proxy_prior, [1 2])
                error(['Value error for proxy prior. Should be 1 or 2.']);
            end
        end        


        function [var_insample_fit] = get_var_insample_fit(self)
            var_insample_fit = self.user_inputs.tab_2_var.insample_fit;
            if ~islogical(var_insample_fit)
                error(['Type error for in-sample fit. Should be boolean.']);
            end
        end 
        

        function [constrained_coefficients] = get_constrained_coefficients(self)
            constrained_coefficients = self.user_inputs.tab_2_var.constrained_coefficients;
            if ~islogical(constrained_coefficients)
                error(['Type error for constrained coefficients. Should be boolean.']);
            end
        end 
        

        function [sums_of_coefficients] = get_sums_of_coefficients(self)
            sums_of_coefficients = self.user_inputs.tab_2_var.sums_of_coefficients;
            if ~islogical(sums_of_coefficients)
                error(['Type error for sums-of-coefficients. Should be boolean.']);
            end
        end          
        

        function [initial_observation] = get_initial_observation(self)
            initial_observation = self.user_inputs.tab_2_var.initial_observation;
            if ~islogical(initial_observation)
                error(['Type error for dummy initial observation. Should be boolean.']);
            end
        end 
        

        function [long_run] = get_long_run(self)
            long_run = self.user_inputs.tab_2_var.long_run;
            if ~islogical(long_run)
                error(['Type error for long-run prior. Should be boolean.']);
            end
        end 
        

        function [stationary] = get_stationary(self)
            stationary = self.user_inputs.tab_2_var.stationary;
            if ~islogical(stationary)
                error(['Type error for stationary prior. Should be boolean.']);
            end
        end 
        

        function [marginal_likelihood] = get_var_marginal_likelihood(self)
            marginal_likelihood = self.user_inputs.tab_2_var.marginal_likelihood;
            if ~islogical(marginal_likelihood)
                error(['Type error for marginal likelihood. Should be boolean.']);
            end
        end
        

        function [hyperparameter_optimization] = get_var_hyperparameter_optimization(self)
            hyperparameter_optimization = self.user_inputs.tab_2_var.hyperparameter_optimization;
            if ~islogical(hyperparameter_optimization)
                error(['Type error for hyperparameter optimization. Should be boolean.']);
            end
        end 
        

        function [coefficients_file] = get_coefficients_file(self)
            coefficients_file = self.user_inputs.tab_2_var.coefficients_file;
            if ~ischar(coefficients_file)
                error(['Type error for constrained coefficients file. Should be char.']);
            end
            coefficients_file = iu.fix_char(coefficients_file);
        end
        

        function [long_run_file] = get_long_run_file(self)
            long_run_file = self.user_inputs.tab_2_var.long_run_file;
            if ~ischar(long_run_file)
                error(['Type error for long-run prior file. Should be char.']);
            end
            long_run_file = iu.fix_char(long_run_file);
        end
        

        function [endogenous, exogenous, dates] = get_var_insample_data(self)
            % check that data path and files are valid
            iu.check_file_path(self.project_path, self.data_file);
            % then load data file
            data = iu.load_data(self.project_path, self.data_file);
            % check that endogenous and exogenous variables are found in data
            iu.check_variables(data, self.data_file, self.endogenous_variables, 'Endogenous variable');
            iu.check_variables(data, self.data_file, self.exogenous_variables, 'Exogenous variable(s)');            
            % check that the start and end dates can be found in the file
            iu.check_dates(data, self.data_file, self.start_date, self.end_date);
            % recover endogenous and exogenous data
            [endogenous] = iu.fetch_data(data, self.data_file, self.start_date, ...
            self.end_date, self.endogenous_variables, 'Endogenous variables');            
            [exogenous] = iu.fetch_data(data, self.data_file, self.start_date, ...
            self.end_date, self.exogenous_variables, 'Exogenous variables');              
            % infer date format, then recover sample dates
            date_format = iu.infer_date_format(self.frequency, ...
                          self.data_file, self.start_date, self.end_date);          
            [dates] = iu.generate_dates(data, date_format, self.frequency, ...
                      self.data_file, self.start_date, self.end_date);    
        end
        

        function [proxys] = get_proxy_data(self)
            % load data only if specified model is proxy-SVAR
            if self.var_type == 7
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.data_file);
                % then load data file
                data = iu.load_data(self.project_path, self.data_file);
                % check that Z variables are found in data
                iu.check_variables(data, self.data_file, self.proxy_variables, 'Proxy variables');
                % check that the start and end dates can be found in the file        
                iu.check_dates(data, self.data_file, self.start_date, self.end_date);
                % recover proxy variables
                [proxys] = iu.fetch_data(data, self.data_file, self.start_date, ...
                         self.end_date, self.proxy_variables, 'Proxy variables');   
            % if model is not proxy-SVAR, return empty list
            else
                proxys = [];
            end
        end 
        

        function [constrained_coefficients_table] = get_constrained_coefficients_table(self)
            % load data only if constrained coefficients is selected
            if self.constrained_coefficients && self.var_type ~= 1
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.coefficients_file);  
                % then load data file
                data = iu.load_data(self.project_path, self.coefficients_file);
                % check data format
                iu.check_coefficients_table(data, self.endogenous_variables, ...
                    self.exogenous_variables, self.lags, self.var_constant, ...
                    self.var_trend, self.var_quadratic_trend, self.coefficients_file);
                constrained_coefficients_table = iu.get_constrained_coefficients_table(data, ...
                    self.endogenous_variables, self.exogenous_variables);
            else
            constrained_coefficients_table = [];
            end
        end


        function [long_run_table] = get_long_run_table(self)
            % load data only if long run prior is selected
            if self.long_run && self.var_type ~= 1
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.long_run_file);
                % then load data file
                data = iu.load_data(self.project_path, self.long_run_file);
                % check data format
                iu.check_long_run_table(data, self.endogenous_variables, self.long_run_file);
                % if format is correct, convert to matrix
                long_run_table = table2array(data);
            else
                long_run_table = [];
            end
        end  


        function [Z_p, Y_p, forecast_dates] = get_var_forecast_data(self)
            % default values for endogenous and exogenous
            Z_p = [];
            Y_p = [];
            % if forecast is selected, recover forecast dates
            if self.forecast || self.conditional_forecast           
                % recover forecast dates, create default values for endogenous and exogenous
                end_date = self.var_dates(end);
                forecast_dates = iu.generate_forecast_dates(end_date, self.forecast_periods, self.frequency);            
            % if forecasts is not selected, return empty data
            else
                forecast_dates = [];
            end
            % if forecast is selected, further recover endogenous and exogenous, if relevant
            if self.forecast && (self.forecast_evaluation || ~iu.is_empty(self.exogenous_variables))
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.forecast_file);
                % then load data file
                data = iu.load_data(self.project_path, self.forecast_file); 
                % if forecast evaluation is selected
                if self.forecast_evaluation
                    % check that exogenous variables are found in data
                    iu.check_variables(data, self.forecast_file, self.endogenous_variables, 'endogenous variables');
                    % load endogenous data
                    Y_p = iu.fetch_forecast_data(data, [], self.endogenous_variables, ...
                    self.forecast_file, self.forecast_evaluation, self.forecast_periods, 'endogenous variable'); 
                end
                % if there are exogenous variables in the model
                if ~iu.is_empty(self.exogenous_variables)
                    % check that exogenous variables are found in data
                    iu.check_variables(data, self.forecast_file, self.exogenous_variables, 'exogenous variables');
                    % load exogenous data 
                    Z_p = iu.fetch_forecast_data(data, [], self.exogenous_variables, ...
                    self.forecast_file, true, self.forecast_periods, 'exogenous variable');
                end
            end
        end


        function [condition_table, shock_table] = get_condition_table(self)
            % if conditional forecast is selected, load data
            if self.conditional_forecast && self.var_type ~= 1
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.conditional_forecast_file);
                % then load data file
                data = iu.load_data(self.project_path, self.conditional_forecast_file);
                % check data format
                iu.check_condition_table(data, self.endogenous_variables, ...
                self.forecast_periods, self.conditional_forecast_file);
                % if format is correct, recover conditions
                [condition_table, shock_table] = iu.get_condition_table(data, self.endogenous_variables);
            % if conditional forecast is not selected, return empty lists
            else
                condition_table = [];
                shock_table = [];
            end
        end
        
        
        function [restriction_table] = get_restriction_table(self)
            % if sign restriction is selected, load data
            if self.structural_identification == 4 && self.var_type ~= 1
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.structural_identification_file);
                % then load data file
                data = iu.load_restriction_data(self.project_path, self.structural_identification_file);
                % get raw sample dates
                raw_dates = iu.get_raw_sample_dates(self.project_path, self.data_file, self.start_date, self.end_date);
                % check data format
                iu.check_restriction_table(data, raw_dates, self.endogenous_variables, self.proxy_variables, ...
                                       self.var_type, self.irf_periods, self.structural_identification_file);
                % if format is correct, recover restrictions
                restriction_table = iu.get_restriction_table(data, raw_dates, self.endogenous_variables, self.proxy_variables);
            % if sign restriction is not selected, return empty list
            else
                restriction_table = [];
            end
        end        
        
    end
    
end
        
        