classdef NowcastingProcessor < handle

    
    properties (GetAccess = public, SetAccess= protected)
        now_model
        now_iterations
        now_burnin
        now_model_credibility
        midas_endogenous_lags
        midas_exogenous_lags
        midas_polynomial_order
        midas_representation
        midas_prior_type
        midas_omega1
        midas_omega2
        midas_upsilon1
        midas_upsilon2
        mfbvar_constant
        mfbvar_trend
        mfbvar_quadratic_trend
        mfbvar_decomposition
        mfbvar_lags
        mfbvar_ar_coefficients
        mfbvar_pi1
        mfbvar_pi2
        mfbvar_pi3
        mfbvar_pi4
        mfbvar_decomposition_file
        dfm_factors
        dfm_loadings_lags
        dfm_factor_lags
        dfm_residual_lags
        dfm_sigma
        dfm_omega
        dfm_delta1
        dfm_pi1
        dfm_pi2
        dfm_pi3
        dfm_omega1
        now_endogenous
        now_exogenous
        now_dates
        now_Z_p
        now_Y_p
        now_forecast_dates
        now_decomposition_table
        now_condition_table
        now_shock_table
        now_restriction_table
    end     
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------     
    
    
    methods (Access = public)
        
        
        function self = NowcastingProcessor()
        end        
        
    
        function nowcasting_inputs(self)

            % recover extension model
            self.now_model = self.get_nowcasting_model();
            % recover iterations
            self.now_iterations = self.get_now_iterations();
            % recover burn-in
            self.now_burnin = self.get_now_burnin();
            % recover credibility level for model estimates
            self.now_model_credibility = self.get_now_model_credibility();
            % recover endogenous lags for midas
            self.midas_endogenous_lags = self.get_midas_endogenous_lags();
            % recover exogenous lags for midas
            self.midas_exogenous_lags = self.get_midas_exogenous_lags();
            % recover polynomial order for midas
            self.midas_polynomial_order = self.get_midas_polynomial_order();
            % recover prior type for midas 
            [self.midas_representation self.midas_prior_type] = self.get_midas_model();
            % recover omega1 for midas 
            self.midas_omega1 = self.get_midas_omega1();
            % recover omega2 for midas 
            self.midas_omega2 = self.get_midas_omega2();
            % recover upsilon1 for midas 
            self.midas_upsilon1 = self.get_midas_upsilon1();
            % recover upsilon2 for midas 
            self.midas_upsilon2 = self.get_midas_upsilon2();
            % recover mfbvar constant
            self.mfbvar_constant = self.get_mfbvar_constant();
            % recover mfbvar trend
            self.mfbvar_trend = self.get_mfbvar_trend();
            % recover mfbvar quadratic trend
            self.mfbvar_quadratic_trend = self.get_mfbvar_quadratic_trend();
            % recover mfbvar decomposition
            self.mfbvar_decomposition = self.get_mfbvar_decomposition();
            % recover mfbvar lags
            self.mfbvar_lags = self.get_mfbvar_lags();    
            % recover AR coefficients
            self.mfbvar_ar_coefficients = self.get_mfbvar_ar_coefficients();
            % recover mfbvar pi1
            self.mfbvar_pi1 = self.get_mfbvar_pi1();
            % recover mfbvar pi2
            self.mfbvar_pi2 = self.get_mfbvar_pi2();
            % recover mfbvar pi3
            self.mfbvar_pi3 = self.get_mfbvar_pi3();
            % recover mfbvar pi4
            self.mfbvar_pi4 = self.get_mfbvar_pi4();
            % get long run prior file
            self.mfbvar_decomposition_file = self.get_mfbvar_decomposition_file();
            % get dfm factors
            self.dfm_factors = self.get_dfm_factors();
            % get dfm loadings lags
            self.dfm_loadings_lags = self.get_dfm_loadings_lags();
            % get dfm factor lags
            self.dfm_factor_lags = self.get_dfm_factor_lags();
            % get dfm residual lags
            self.dfm_residual_lags = self.get_dfm_residual_lags();
            % get dfm sigma
            self.dfm_sigma = self.get_dfm_sigma();
            % get dfm omega
            self.dfm_omega = self.get_dfm_omega();
            % get dfm delta1
            self.dfm_delta1 = self.get_dfm_delta1();            
            % get dfm pi1
            self.dfm_pi1 = self.get_dfm_pi1();
            % get dfm pi2
            self.dfm_pi2 = self.get_dfm_pi2();
            % get dfm pi3
            self.dfm_pi3 = self.get_dfm_pi3();
            % get dfm omega1
            self.dfm_omega1 = self.get_dfm_omega1();      
        end
        

        function nowcasting_data(self)
            % print loading message
            if self.progress_bar
                cu.print_message_to_complete('Data loading:');
            end
            % recover in-sample endogenous and exogenous
            [self.now_endogenous, self.now_exogenous, self.now_dates] = self.get_now_insample_data();
            % recover forecast data
            [self.now_Z_p, self.now_Y_p, self.now_forecast_dates] = self.get_now_forecast_data();
            % recover decomposition data
            [self.now_decomposition_table] = self.get_now_decomposition_table();            
            % recover conditional forecast data
            [self.now_condition_table, self.now_shock_table] = self.get_now_condition_table();
            % recover sign restrictions data
            [self.now_restriction_table] = self.get_now_restriction_table();
            % print loading done message
            if self.progress_bar
                cu.print_message('  —  done');
            end
        end


        function make_nowcasting_information(self)
            % get sample dates
            self.results_information.dates = self.now_dates;
            % get forecast dates
            self.results_information.forecast_dates = self.now_forecast_dates;
            self.results_information.conditional_forecast_dates = self.now_forecast_dates;
            % get decomposition file
            self.results_information.decomposition_file = self.mfbvar_decomposition_file;
        end


        function make_nowcasting_graphics_information(self)
            % get sample dates
            self.graphics_information.dates = self.now_dates;
            % get forecast dates
            self.graphics_information.forecast_dates = self.now_forecast_dates;
            self.graphics_information.conditional_forecast_dates = self.now_forecast_dates;
            % get actual data for forecast evaluation, if available
            self.graphics_information.Y_p = self.now_Y_p;
        end


    end
        

    %---------------------------------------------------
    % Methods (Access = private)
    %--------------------------------------------------- 
        
  
    methods (Access = protected, Hidden = true) 
    

        function [model] = get_nowcasting_model(self)
            model = self.user_inputs.tab_2_now.model;
            if ~ismember(model, [1 2 3])
                error(['Value error for nowcasting model. Should be 1, 2 or 3.']);
            end
        end


        function [iterations] = get_now_iterations(self)
            iterations = self.user_inputs.tab_2_now.iterations;
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


        function [burnin] = get_now_burnin(self)
            burnin = self.user_inputs.tab_2_now.burnin;
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


        function [model_credibility] =  get_now_model_credibility(self)
            model_credibility = self.user_inputs.tab_2_now.model_credibility;
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


        function [lags] = get_midas_endogenous_lags(self)
            lags = self.user_inputs.tab_2_now.midas_endogenous_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for MIDAS endogenous lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for MIDAS endogenous lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags < 0
                error(['Value error for MIDAS endogenous lags. Should be positive integer.']);
            end            
        end 


        function [lags] = get_midas_exogenous_lags(self)
            lags = self.user_inputs.tab_2_now.midas_exogenous_lags;
            if ~ismember(class(lags), ["char" "double"])
                error(['Type error for MIDAS exogenous lags. Should be integer or list of integers.']);
            end
            if ischar(lags)
                lags = iu.char_to_array(lags);
                if any(isnan(str2double(lags)))
                    error(['Type error for exogenous lags. All elements should be integers.']);
                else
                    lags = str2double(lags)';
                end
            end
            if numel(lags) ~= numel(self.exogenous_variables) && numel(lags) ~= 1
                error(['Dimension error for exogenous lags. Dimension of exogenous lags and exogenous variables don''t match.']);
            end
            if any(isnan(lags)) || any(isinf(lags))
                error(['Type error for exogenous lags. All elements should be integers.']);
            end
        end
            
        
        function [order] = get_midas_polynomial_order(self)
            order = self.user_inputs.tab_2_now.midas_polynomial_order;
            if ~(ischar(order) || iu.is_integer(order))
                error(['Type error for MIDAS polynomial order. Should be integer.']);
            end
            if ~iu.is_empty(order) && ischar(order)
                if iu.is_digit(order)
                    order = str2double(order);
                else
                    error(['Type error for MIDAS polynomial order. Should be positive integer.']);
                end
            end  
            if iu.is_integer(order) && order <= 0
                error(['Value error for MIDAS polynomial order. Should be positive integer.']);
            end            
        end 


        function [representation prior_type] = get_midas_model(self)
            prior = self.user_inputs.tab_2_now.midas_model;
            if ~ismember(prior, [1 2 3 4 5 6 7 8 9])
                error(['Value error for MIDAS prior. Should be integer between 1 and 9.']);
            end
            if prior == 1
                representation = 'unrestricted';
                prior_type = 'minnesota';
            elseif prior == 2
                representation = 'unrestricted';                
                prior_type = 'horseshoe';
            elseif prior == 3
                representation = 'unrestricted';                
                prior_type = 'lasso';
            elseif prior == 4
                representation = 'almon';                
                prior_type = 'minnesota';
            elseif prior == 5
                representation = 'almon';                
                prior_type = 'horseshoe';
            elseif prior == 6
                representation = 'almon';                
                prior_type = 'lasso';
            elseif prior == 7
                representation = 'fourier';                
                prior_type = 'minnesota';
            elseif prior == 8
                representation = 'fourier';                
                prior_type = 'horseshoe';
            elseif prior == 9
                representation = 'fourier';                
                prior_type = 'lasso';
            end
        end


        function [omega1] = get_midas_omega1(self)
            omega1 = self.user_inputs.tab_2_now.midas_omega1;
            if ~ismember(class(omega1), ["char" "double"])
                error(['Type error for MIDAS omega1. Should be float or integer.']);
            end
            if ischar(omega1)
                if isnan(str2double(omega1))
                    error(['Type error for MIDAS omega1. Should be float or integer.']);
                else
                    omega1 = str2double(omega1);
                end
            end
            if omega1 <= 0
                error(['Value error for MIDAS omega1. Should be strictly positive.']);
            end
        end  


        function [omega2] = get_midas_omega2(self)
            omega2 = self.user_inputs.tab_2_now.midas_omega2;
            if ~ismember(class(omega2), ["char" "double"])
                error(['Type error for MIDAS omega2. Should be float or integer.']);
            end
            if ischar(omega2)
                if isnan(str2double(omega2))
                    error(['Type error for MIDAS omega2. Should be float or integer.']);
                else
                    omega2 = str2double(omega2);
                end
            end
            if omega2 <= 0
                error(['Value error for MIDAS omega2. Should be strictly positive.']);
            end
        end 


        function [upsilon1] = get_midas_upsilon1(self)
            upsilon1 = self.user_inputs.tab_2_now.midas_upsilon1;
            if ~ismember(class(upsilon1), ["char" "double"])
                error(['Type error for MIDAS upsilon1. Should be float or integer.']);
            end
            if ischar(upsilon1)
                if isnan(str2double(upsilon1))
                    error(['Type error for MIDAS upsilon1. Should be float or integer.']);
                else
                    upsilon1 = str2double(upsilon1);
                end
            end
            if upsilon1 <= 0
                error(['Value error for MIDAS upsilon1. Should be strictly positive.']);
            end
        end  


        function [upsilon2] = get_midas_upsilon2(self)
            upsilon2 = self.user_inputs.tab_2_now.midas_upsilon2;
            if ~ismember(class(upsilon2), ["char" "double"])
                error(['Type error for MIDAS upsilon2. Should be float or integer.']);
            end
            if ischar(upsilon2)
                if isnan(str2double(upsilon2))
                    error(['Type error for MIDAS upsilon2. Should be float or integer.']);
                else
                    upsilon2 = str2double(upsilon2);
                end
            end
            if upsilon2 <= 0
                error(['Value error for MIDAS upsilon2. Should be strictly positive.']);
            end
        end 


        function [constant] = get_mfbvar_constant(self)
            constant = self.user_inputs.tab_2_now.mfbvar_constant;
            if ~islogical(constant)
                error(['Type error for constant. Should be boolean.']);
            end
        end         
        

        function [trend] = get_mfbvar_trend(self)
            trend = self.user_inputs.tab_2_now.mfbvar_trend;
            if ~islogical(trend)
                error(['Type error for trend. Should be boolean.']);
            end
        end         
        

        function [quadratic_trend] = get_mfbvar_quadratic_trend(self)
            quadratic_trend = self.user_inputs.tab_2_now.mfbvar_quadratic_trend;
            if ~islogical(quadratic_trend)
                error(['Type error for quadratic trend. Should be boolean.']);
            end
        end 


        function [decomposition] = get_mfbvar_decomposition(self)
            decomposition = self.user_inputs.tab_2_now.mfbvar_decomposition;
            if ~islogical(decomposition)
                error(['Type error for decomposition. Should be boolean.']);
            end
        end 


        function [lags] = get_mfbvar_lags(self)
            lags = self.user_inputs.tab_2_now.mfbvar_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for MF-BVAR lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for MF-BVAR lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags <= 0
                error(['Value error for MF-BVAR lags. Should be positive integer.']);
            end            
        end          
        

        function [ar_coefficients] = get_mfbvar_ar_coefficients(self)
            ar_coefficients = self.user_inputs.tab_2_now.mfbvar_ar_coefficients;
            if ~ismember(class(ar_coefficients), ["char" "double"])
                error(['Type error for MF-BVAR AR coefficients. Should be float or scalar array.']);
            end
            if ischar(ar_coefficients)
                ar_coefficients = iu.char_to_array(ar_coefficients);
                if any(isnan(str2double(ar_coefficients)))
                    error(['Type error for MF-BVAR AR coefficients. All elements should be scalars.']);
                else
                    ar_coefficients = str2double(ar_coefficients)';
                end
            end
            if numel(ar_coefficients) ~= numel(self.endogenous_variables) && numel(ar_coefficients) ~= 1
                error(['Dimension error for MF-BVAR AR coefficients. Dimension of AR coefficients and endogenous don''t match.']);
            end
            if any(isnan(ar_coefficients)) || any(isinf(ar_coefficients))
                error(['Type error for MF-BVAR AR coefficients. All elements should be scalars.']);
            end
        end
        

        function [pi1] = get_mfbvar_pi1(self)
            pi1 = self.user_inputs.tab_2_now.mfbvar_pi1;
            if ~ismember(class(pi1), ["char" "double"])
                error(['Type error for MF-BVAR pi1. Should be float or integer.']);
            end
            if ischar(pi1)
                if isnan(str2double(pi1))
                    error(['Type error for MF-BVAR pi1. Should be float or integer.']);
                else
                    pi1 = str2double(pi1);
                end
            end
            if pi1 <= 0
                error(['Value error for MF-BVAR pi1. Should be strictly positive.']);
            end
        end   
        

        function [pi2] = get_mfbvar_pi2(self)
            pi2 = self.user_inputs.tab_2_now.mfbvar_pi2;
            if ~ismember(class(pi2), ["char" "double"])
                error(['Type error for MF-BVAR pi2. Should be float or integer.']);
            end
            if ischar(pi2)
                if isnan(str2double(pi2))
                    error(['Type error for MF-BVAR pi2. Should be float or integer.']);
                else
                    pi2 = str2double(pi2);
                end
            end
            if pi2 <= 0
                error(['Value error for MF-BVAR pi2. Should be strictly positive.']);
            end
        end   
        

        function [pi3] = get_mfbvar_pi3(self)
            pi3 = self.user_inputs.tab_2_now.mfbvar_pi3;
            if ~ismember(class(pi3), ["char" "double"])
                error(['Type error for MF-BVAR pi3. Should be float or integer.']);
            end
            if ischar(pi3)
                if isnan(str2double(pi3))
                    error(['Type error for MF-BVAR pi3. Should be float or integer.']);
                else
                    pi3 = str2double(pi3);
                end
            end
            if pi3 <= 0
                error(['Value error for MF-BVAR pi3. Should be strictly positive.']);
            end
        end          
        

        function [pi4] = get_mfbvar_pi4(self)
            pi4 = self.user_inputs.tab_2_now.mfbvar_pi4;
            if ~ismember(class(pi4), ["char" "double"])
                error(['Type error for MF-BVAR pi4. Should be float or integer.']);
            end
            if ischar(pi4)
                if isnan(str2double(pi4))
                    error(['Type error for MF-BVAR pi4. Should be float or integer.']);
                else
                    pi4 = str2double(pi4);
                end
            end
            if pi4 <= 0
                error(['Value error for MF-BVAR pi4. Should be strictly positive.']);
            end
        end         


        function [decomposition_file] = get_mfbvar_decomposition_file(self)
            decomposition_file = self.user_inputs.tab_2_now.mfbvar_decomposition_file;
            if ~ischar(decomposition_file)
                error(['Type error for MF-BVAR decomposition file. Should be char.']);
            end
            decomposition_file = iu.fix_char(decomposition_file);
        end


        function [factors] = get_dfm_factors(self)
            factors = self.user_inputs.tab_2_now.dfm_factors;
            if ~(ischar(factors) || iu.is_integer(factors))
                error(['Type error for DFM factors. Should be integer.']);
            end
            if ~iu.is_empty(factors) && ischar(factors)
                if iu.is_digit(factors)
                    factors = str2double(factors);
                else
                    error(['Type error for DFM factors. Should be positive integer.']);
                end
            end  
            if iu.is_integer(factors) && factors <= 0
                error(['Value error for DFM factors. Should be positive integer.']);
            end            
        end 


        function [lags] = get_dfm_loadings_lags(self)
            lags = self.user_inputs.tab_2_now.dfm_loadings_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for DFM loadings lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for DFM loadings lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags < 0
                error(['Value error for DFM loadings lags. Should be positive integer.']);
            end            
        end 


        function [lags] = get_dfm_factor_lags(self)
            lags = self.user_inputs.tab_2_now.dfm_factor_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for DFM factor lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for DFM factor lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags <= 0
                error(['Value error for DFM factor lags. Should be positive integer.']);
            end            
        end 


        function [lags] = get_dfm_residual_lags(self)
            lags = self.user_inputs.tab_2_now.dfm_residual_lags;
            if ~(ischar(lags) || iu.is_integer(lags))
                error(['Type error for DFM residual lags. Should be integer.']);
            end
            if ~iu.is_empty(lags) && ischar(lags)
                if iu.is_digit(lags)
                    lags = str2double(lags);
                else
                    error(['Type error for DFM residual lags. Should be positive integer.']);
                end
            end  
            if iu.is_integer(lags) && lags < 0
                error(['Value error for DFM residual lags. Should be positive integer.']);
            end            
        end 


        function [sigma] = get_dfm_sigma(self)
            sigma = self.user_inputs.tab_2_now.dfm_sigma;
            if ~ismember(class(sigma), ["char" "double"])
                error(['Type error for DFM sigma. Should be float or integer.']);
            end
            if ischar(sigma)
                if isnan(str2double(sigma))
                    error(['Type error for DFM sigma. Should be float or integer.']);
                else
                    sigma = str2double(sigma);
                end
            end
            if sigma <= 0
                error(['Value error for DFM sigma. Should be strictly positive.']);
            end
        end  


        function [omega] = get_dfm_omega(self)
            omega = self.user_inputs.tab_2_now.dfm_omega;
            if ~ismember(class(omega), ["char" "double"])
                error(['Type error for DFM omega. Should be float or integer.']);
            end
            if ischar(omega)
                if isnan(str2double(omega))
                    error(['Type error for DFM omega. Should be float or integer.']);
                else
                    omega = str2double(omega);
                end
            end
            if omega <= 0
                error(['Value error for DFM omega. Should be strictly positive.']);
            end
        end  


        function [delta1] = get_dfm_delta1(self)
            delta1 = self.user_inputs.tab_2_now.dfm_delta1;
            if ~ismember(class(delta1), ["char" "double"])
                error(['Type error for DFM delta1. Should be float or integer.']);
            end
            if ischar(delta1)
                if isnan(str2double(delta1))
                    error(['Type error for DFM delta1. Should be float or integer.']);
                else
                    delta1 = str2double(delta1);
                end
            end
            if delta1 <= 0
                error(['Value error for DFM delta1. Should be strictly positive.']);
            end
        end  


        function [pi1] = get_dfm_pi1(self)
            pi1 = self.user_inputs.tab_2_now.dfm_pi1;
            if ~ismember(class(pi1), ["char" "double"])
                error(['Type error for DFM pi1. Should be float or integer.']);
            end
            if ischar(pi1)
                if isnan(str2double(pi1))
                    error(['Type error for DFM pi1. Should be float or integer.']);
                else
                    pi1 = str2double(pi1);
                end
            end
            if pi1 <= 0
                error(['Value error for DFM pi1. Should be strictly positive.']);
            end
        end   


        function [pi2] = get_dfm_pi2(self)
            pi2 = self.user_inputs.tab_2_now.dfm_pi2;
            if ~ismember(class(pi2), ["char" "double"])
                error(['Type error for DFM pi2. Should be float or integer.']);
            end
            if ischar(pi2)
                if isnan(str2double(pi2))
                    error(['Type error for DFM pi2. Should be float or integer.']);
                else
                    pi2 = str2double(pi2);
                end
            end
            if pi2 <= 0
                error(['Value error for DFM pi2. Should be strictly positive.']);
            end
        end 


        function [pi3] = get_dfm_pi3(self)
            pi3 = self.user_inputs.tab_2_now.dfm_pi3;
            if ~ismember(class(pi3), ["char" "double"])
                error(['Type error for DFM pi3. Should be float or integer.']);
            end
            if ischar(pi3)
                if isnan(str2double(pi3))
                    error(['Type error for DFM pi3. Should be float or integer.']);
                else
                    pi3 = str2double(pi3);
                end
            end
            if pi3 <= 0
                error(['Value error for DFM pi3. Should be strictly positive.']);
            end
        end 


        function [omega1] = get_dfm_omega1(self)
            omega1 = self.user_inputs.tab_2_now.dfm_omega1;
            if ~ismember(class(omega1), ["char" "double"])
                error(['Type error for DFM omega1. Should be float or integer.']);
            end
            if ischar(omega1)
                if isnan(str2double(omega1))
                    error(['Type error for DFM omega1. Should be float or integer.']);
                else
                    omega1 = str2double(omega1);
                end
            end
            if omega1 <= 0
                error(['Value error for DFM omega1. Should be strictly positive.']);
            end
        end 


        function [endogenous, exogenous, dates] = get_now_insample_data(self)
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
            [endogenous] = iu.fetch_nan_data(data, self.data_file, self.start_date, ...
            self.end_date, self.endogenous_variables, 'Endogenous variables');            
            [exogenous] = iu.fetch_nan_data(data, self.data_file, self.start_date, ...
            self.end_date, self.exogenous_variables, 'Exogenous variables');              
            % infer date format, then recover sample dates
            date_format = iu.infer_date_format(self.frequency, ...
                          self.data_file, self.start_date, self.end_date);          
            [dates] = iu.generate_dates(data, date_format, self.frequency, ...
                      self.data_file, self.start_date, self.end_date);   
            if self.now_model == 3
                start_date = find(strcmp(data.Properties.RowNames, self.start_date));
                end_date = find(strcmp(data.Properties.RowNames, self.end_date));
                sample_data = data(start_date:end_date,:);
                dates = datetime(sample_data.Properties.RowNames(~isnan ...
                        (sample_data.(convertStringsToChars(self.endogenous_variables)))));
            end
        end


        function [Z_p, Y_p, forecast_dates] = get_now_forecast_data(self)
            % default values for endogenous and exogenous
            Z_p = [];
            Y_p = [];
            % if forecast is selected, recover forecast dates
            if self.forecast || self.conditional_forecast           
                % recover forecast dates, create default values for endogenous and exogenous
                end_date = self.now_dates(end);
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
                % if the model if MFBVAR and there are exogenous variables in the model
                if self.now_model == 1 && ~iu.is_empty(self.exogenous_variables)
                    % load in-sample data to fill missing variables
                    in_sample_data = iu.load_data(self.project_path, self.data_file); 
                    in_sample_data = table2array(in_sample_data(end-1:end,self.exogenous_variables));
                    % load exogenous data 
                    Z_p = iu.fetch_forecast_data(data, in_sample_data, self.exogenous_variables, ...
                    self.forecast_file, true, self.forecast_periods, 'exogenous variable');
                end
            end
        end


        function [decomposition_table] = get_now_decomposition_table(self)
            % if model is MFBVAR and decomposition is selected, load data
            if self.now_model == 1 && self.mfbvar_decomposition
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.mfbvar_decomposition_file);
                % then load data file
                data = iu.load_data(self.project_path, self.mfbvar_decomposition_file);
                % check data format
                iu.check_decomposition_table(data, self.endogenous_variables, self.mfbvar_decomposition_file);
                % turn to matrix
                decomposition_table = table2array(data);
            else
                decomposition_table = [];
            end
        end


        function [condition_table, shock_table] = get_now_condition_table(self)
            % if model is MFBVAR and conditional forecast is selected, load data
            if self.now_model == 1 && self.conditional_forecast
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


        function [restriction_table] = get_now_restriction_table(self)
            % if model is MFBVAR and sign restriction is selected, load data
            if self.now_model == 1 && self.structural_identification == 4
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
        
        