classdef VecVarmaProcessor < handle

    
    properties (GetAccess = public, SetAccess= protected)

        ext_model
        ext_iterations
        ext_burnin
        ext_model_credibility
        ext_constant
        ext_trend
        ext_quadratic_trend
        vec_lags
        vec_pi1
        vec_pi2
        vec_pi3
        vec_pi4
        vec_prior_type
        error_correction_type
        max_cointegration_rank
        varma_lags
        varma_ar_coefficients
        varma_pi1
        varma_pi2
        varma_pi3
        varma_pi4
        residual_lags
        lambda1
        lambda2
        lambda3
        ext_endogenous
        ext_exogenous
        ext_dates
        ext_Z_p
        ext_Y_p
        ext_forecast_dates
        ext_condition_table
        ext_shock_table
        ext_restriction_table
    end     
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------     
    
    
    methods (Access = public)
        
        
        function self = VecVarmaProcessor()
        end        
        
    
        function vec_varma_inputs(self)

            % recover extension model
            self.ext_model = self.get_extension_model();
            % recover iterations
            self.ext_iterations = self.get_ext_iterations();
            % recover burn-in
            self.ext_burnin = self.get_ext_burnin();
            % recover credibility level for model estimates
            self.ext_model_credibility = self.get_ext_model_credibility();
            % recover constant
            self.ext_constant = self.get_ext_constant(); 
            % recover trend
            self.ext_trend = self.get_ext_trend(); 
            % recover quadratic trend
            self.ext_quadratic_trend = self.get_ext_quadratic_trend(); 
            % recover lags
            self.vec_lags = self.get_vec_lags();         
            % recover vec pi1
            self.vec_pi1 = self.get_vec_pi1();  
            % recover vec pi2
            self.vec_pi2 = self.get_vec_pi2(); 
            % recover vec pi3
            self.vec_pi3 = self.get_vec_pi3(); 
            % recover vec pi4
            self.vec_pi4 = self.get_vec_pi4(); 
            % recover prior type
            self.vec_prior_type = self.get_prior_type();
            % recover error correction type
            self.error_correction_type = self.get_error_correction_type();
            % recover max cointegration rank
            self.max_cointegration_rank = self.get_max_cointegration_rank();
            % recover VARMA lags
            self.varma_lags = self.get_varma_lags();
            % recover AR coefficients
            self.varma_ar_coefficients = self.get_varma_ar_coefficients(); 
            % recover varma pi1
            self.varma_pi1 = self.get_varma_pi1();  
            % recover varma pi2
            self.varma_pi2 = self.get_varma_pi2(); 
            % recover varma pi3
            self.varma_pi3 = self.get_varma_pi3(); 
            % recover varma pi4
            self.varma_pi4 = self.get_varma_pi4(); 
            % recover residual lags
            self.residual_lags = self.get_residual_lags();  
            % recover lambda1
            self.lambda1 = self.get_lambda1();
            % recover lambda2
            self.lambda2 = self.get_lambda2();
            % recover lambda3
            self.lambda3 = self.get_lambda3();       
        end
        

        function vec_varma_data(self)
            % print loading message
            if self.progress_bar
                cu.print_message_to_complete('Data loading:');
            end
            % recover in-sample endogenous and exogenous
            [self.ext_endogenous, self.ext_exogenous, self.ext_dates] = self.get_ext_insample_data();
            % recover forecast data
            [self.ext_Z_p, self.ext_Y_p, self.ext_forecast_dates] = self.get_ext_forecast_data();
            % recover conditional forecast data
            [self.ext_condition_table, self.ext_shock_table] = self.get_ext_condition_table();
            % recover sign restrictions data
            [self.ext_restriction_table] = self.get_ext_restriction_table();
            if self.progress_bar
                cu.print_message('  —  done');
            end
        end


        function make_vec_varma_information(self)
            % get sample dates
            self.results_information.dates = self.ext_dates;
            % get forecast dates
            self.results_information.forecast_dates = self.ext_forecast_dates;
            self.results_information.conditional_forecast_dates = self.ext_forecast_dates;
        end


        function make_vec_varma_graphics_information(self)
            % get sample dates
            self.graphics_information.dates = self.ext_dates;
            % get forecast dates
            self.graphics_information.forecast_dates = self.ext_forecast_dates;
            self.graphics_information.conditional_forecast_dates = self.ext_forecast_dates;
            % get actual data for forecast evaluation, if available
            self.graphics_information.Y_p = self.ext_Y_p;
        end


    end
        

    %---------------------------------------------------
    % Methods (Access = private)
    %--------------------------------------------------- 
        
  
    methods (Access = protected, Hidden = true) 
    

        function [model] = get_extension_model(self)
            model = self.user_inputs.tab_2_ext.model;
            if ~ismember(model, [1 2])
                error(['Value error for VAR extension type. Should be 1 or 2.']);
            end
        end


        function [iterations] = get_ext_iterations(self)
            iterations = self.user_inputs.tab_2_ext.iterations;
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


        function [burnin] = get_ext_burnin(self)
            burnin = self.user_inputs.tab_2_ext.burnin;
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


        function [model_credibility] =  get_ext_model_credibility(self)
            model_credibility = self.user_inputs.tab_2_ext.model_credibility;
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


        function [constant] = get_ext_constant(self)
            constant = self.user_inputs.tab_2_ext.constant;
            if ~islogical(constant)
                error(['Type error for constant. Should be boolean.']);
            end
        end         
        

        function [trend] = get_ext_trend(self)
            trend = self.user_inputs.tab_2_ext.trend;
            if ~islogical(trend)
                error(['Type error for trend. Should be boolean.']);
            end
        end         
        

        function [quadratic_trend] = get_ext_quadratic_trend(self)
            quadratic_trend = self.user_inputs.tab_2_ext.quadratic_trend;
            if ~islogical(quadratic_trend)
                error(['Type error for quadratic trend. Should be boolean.']);
            end
        end 
        

        function [lags] = get_vec_lags(self)
            lags = self.user_inputs.tab_2_ext.vec_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for VEC lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for VEC lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags <= 0
                error(['Value error for VEC lags. Should be positive integer.']);
            end            
        end        


        function [pi1] = get_vec_pi1(self)
            pi1 = self.user_inputs.tab_2_ext.vec_pi1;
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
        

        function [pi2] = get_vec_pi2(self)
            pi2 = self.user_inputs.tab_2_ext.vec_pi2;
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
        

        function [pi3] = get_vec_pi3(self)
            pi3 = self.user_inputs.tab_2_ext.vec_pi3;
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
        

        function [pi4] = get_vec_pi4(self)
            pi4 = self.user_inputs.tab_2_ext.vec_pi4;
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
        

        function [prior] = get_prior_type(self)
            prior = self.user_inputs.tab_2_ext.prior_type;
            if ~ismember(prior, [1 2 3])
                error(['Value error for VEC prior type. Should be 1, 2 or 3.']);
            end
        end

        
        function [error_correction_type] = get_error_correction_type(self)
            error_correction_type = self.user_inputs.tab_2_ext.error_correction_type;
            if ~ismember(error_correction_type, [1 2])
                error(['Value error for VEC error correction type. Should be 1 or 2.']);
            end
        end


        function [rank] = get_max_cointegration_rank(self)
            rank = self.user_inputs.tab_2_ext.max_cointegration_rank;
            if ~(ischar(rank) || iu.is_integer(rank))
                error(['Type error for max cointegration rank. Should be integer.']);
            end
            if ~iu.is_empty(rank) && ischar(rank)
                if iu.is_digit(rank)
                    rank = str2double(rank);
                else
                    error(['Type error for max cointegration rank. Should be positive integer.']);
                end
            end  
            if iu.is_integer(rank) && rank <= 0
                error(['Value error for max cointegration rank. Should be positive integer.']);
            end            
        end


        function [lags] = get_varma_lags(self)
            lags = self.user_inputs.tab_2_ext.varma_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for VARMA lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for VARMA lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags <= 0
                error(['Value error for VARMA lags. Should be positive integer.']);
            end            
        end


        function [ar_coefficients] = get_varma_ar_coefficients(self)
            ar_coefficients = self.user_inputs.tab_2_ext.ar_coefficients;
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


        function [pi1] = get_varma_pi1(self)
            pi1 = self.user_inputs.tab_2_ext.varma_pi1;
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
        

        function [pi2] = get_varma_pi2(self)
            pi2 = self.user_inputs.tab_2_ext.varma_pi2;
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
        

        function [pi3] = get_varma_pi3(self)
            pi3 = self.user_inputs.tab_2_ext.varma_pi3;
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
        

        function [pi4] = get_varma_pi4(self)
            pi4 = self.user_inputs.tab_2_ext.varma_pi4;
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
        

        function [lags] = get_residual_lags(self)
            lags = self.user_inputs.tab_2_ext.residual_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for residual lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for residual lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags <= 0
                error(['Value error for residual lags. Should be positive integer.']);
            end            
        end


        function [lambda1] = get_lambda1(self)
            lambda1 = self.user_inputs.tab_2_ext.lambda1;
            if ~ismember(class(lambda1), ["char" "double"])
                error(['Type error for lambda1. Should be float or integer.']);
            end
            if ischar(lambda1)
                if isnan(str2double(lambda1))
                    error(['Type error for lambda1. Should be float or integer.']);
                else
                    lambda1 = str2double(lambda1);
                end
            end
            if lambda1 <= 0
                error(['Value error for lambda1. Should be strictly positive.']);
            end
        end   
        

        function [lambda2] = get_lambda2(self)
            lambda2 = self.user_inputs.tab_2_ext.lambda2;
            if ~ismember(class(lambda2), ["char" "double"])
                error(['Type error for lambda2. Should be float or integer.']);
            end
            if ischar(lambda2)
                if isnan(str2double(lambda2))
                    error(['Type error for lambda2. Should be float or integer.']);
                else
                    lambda2 = str2double(lambda2);
                end
            end
            if lambda2 <= 0
                error(['Value error for lambda2. Should be strictly positive.']);
            end
        end  
        

        function [lambda3] = get_lambda3(self)
            lambda3 = self.user_inputs.tab_2_ext.lambda3;
            if ~ismember(class(lambda3), ["char" "double"])
                error(['Type error for lambda3. Should be float or integer.']);
            end
            if ischar(lambda3)
                if isnan(str2double(lambda3))
                    error(['Type error for lambda3. Should be float or integer.']);
                else
                    lambda3 = str2double(lambda3);
                end
            end
            if lambda3 <= 0
                error(['Value error for lambda3. Should be strictly positive.']);
            end
        end  


        function [endogenous, exogenous, dates] = get_ext_insample_data(self)
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


        function [Z_p, Y_p, forecast_dates] = get_ext_forecast_data(self)
            % default values for endogenous and exogenous
            Z_p = [];
            Y_p = [];
            % if forecast is selected, recover forecast dates
            if self.forecast || self.conditional_forecast           
                % recover forecast dates, create default values for endogenous and exogenous
                end_date = self.ext_dates(end);
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


        function [condition_table, shock_table] = get_ext_condition_table(self)
            % if conditional forecast is selected, load data
            if self.conditional_forecast
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
        
        
        function [restriction_table] = get_ext_restriction_table(self)
            % if sign restriction is selected, load data
            if self.structural_identification == 4
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.structural_identification_file);
                % then load data file
                data = iu.load_restriction_data(self.project_path, self.structural_identification_file);
                % get raw sample dates
                raw_dates = iu.get_raw_sample_dates(self.project_path, self.data_file, self.start_date, self.end_date);
                % check data format
                iu.check_restriction_table(data, raw_dates, self.endogenous_variables, [], ...
                                           2, self.irf_periods, self.structural_identification_file);
                % if format is correct, recover restrictions
                restriction_table = iu.get_restriction_table(data, raw_dates, self.endogenous_variables, []);
            % if sign restriction is not selected, return empty list
            else
                restriction_table = [];
            end
        end      
   
        
    end
    
end
        
        