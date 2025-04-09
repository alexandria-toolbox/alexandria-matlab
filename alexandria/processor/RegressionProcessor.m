classdef RegressionProcessor < handle

    
    properties (GetAccess = public, SetAccess= protected)
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
        endogenous
        exogenous
        dates
        Z
        X_p
        y_p
        Z_p
        forecast_dates
    end     
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------     
    
    
    methods (Access = public)
        
        
        function self = RegressionProcessor()
        end        
        
    
        function regression_inputs(self)
            % recover regression type
            self.regression_type = self.get_regression_type();
            % recover iterations
            self.iterations = self.get_iterations();
            % recover burn-in
            self.burnin = self.get_burnin();
            % recover credibility level for model estimates
            self.model_credibility = self.get_model_credibility();
            % recover b
            self.b = self.get_b();   
            % recover V
            self.V = self.get_V();
            % recover alpha
            self.alpha = self.get_alpha();             
            % recover alpha
            self.delta = self.get_delta();
            % recover g
            self.g = self.get_g();          
            % recover Q
            self.Q = self.get_Q();
            % recover tau
            self.tau = self.get_tau();
            % recover thinning
            self.thinning = self.get_thinning();     
            % recover thinning frequency
            self.thinning_frequency = self.get_thinning_frequency();            
            % recover bame of file for Z regressors
            self.Z_variables = self.get_Z_variables();
            % recover q
            self.q = self.get_q(); 
            % recover p
            self.p = self.get_p(); 
            % recover H
            self.H = self.get_H();             
            % recover constant
            self.constant = self.get_constant();            
            % recover b (constant)
            self.b_constant = self.get_b_constant();
            % recover V (constant)
            self.V_constant = self.get_V_constant();         
            % recover trend
            self.trend = self.get_trend();            
            % recover b (trend)
            self.b_trend = self.get_b_trend();
            % recover V (constant)
            self.V_trend = self.get_V_trend();              
            % recover quadratic trend
            self.quadratic_trend = self.get_quadratic_trend();
            % recover b (quadratic trend)
            self.b_quadratic_trend = self.get_b_quadratic_trend();
            % recover V (quadratic trend)
            self.V_quadratic_trend = self.get_V_quadratic_trend();           
            % recover in-sample fit
            self.insample_fit = self.get_insample_fit(); 
            % recover marginal likelihood
            self.marginal_likelihood = self.get_marginal_likelihood(); 
            % recover hyperparameter optimization
            self.hyperparameter_optimization = self.get_hyperparameter_optimization();             
            % recover optimization type
            self.optimization_type = self.get_optimization_type();        
        end
        
        
        function regression_data(self)
            % print loading message
            if self.progress_bar
                cu.print_message('Data loading:');
            end
            % recover in-sample endogenous and exogenous
            [self.endogenous, self.exogenous, self.dates] = self.get_insample_data();
            % recover heteroscedastic data
            self.Z = self.get_heteroscedastic_data();
            % recover forecast data
            [self.X_p, self.y_p, self.Z_p, self.forecast_dates] = self.get_forecast_data();
            % print loading completion message
            if self.progress_bar            
                cu.print_message(['  — / —    [' repmat('=',1,33) ']  —  done']);
            end
        end

    
        function make_regression_information(self)
            % get sample dates
            self.results_information.dates = self.dates;       
            % get heteroscedastic variables
            if self.regression_type == 5
                self.results_information.heteroscedastic_variables = self.Z_variables;
            else
                self.results_information.heteroscedastic_variables = [];
            end
            % get regression option: in-sample fit
            self.results_information.insample_fit = self.insample_fit;
            % get regression option: marginal likelihood
            self.results_information.marginal_likelihood = self.marginal_likelihood;
            % get regression option: hyperparameter optimization
            self.results_information.hyperparameter_optimization = self.hyperparameter_optimization;
            % get regression option: optimization type
            if self.optimization_type == 1
                self.results_information.optimization_type = 'simple';
            elseif self.optimization_type == 2
                self.results_information.optimization_type = 'full';
            end
        end

    
        function make_regression_graphics_information(self)
            % get sample dates
            self.graphics_information.dates = self.dates;
            % get forecast dates
            self.graphics_information.forecast_dates = self.forecast_dates;
            % get actual data for forecast evaluation, if available
            self.graphics_information.y_p = self.y_p;
        end


    end  
    
    
    %---------------------------------------------------
    % Methods (Access = private)
    %--------------------------------------------------- 
        
  
    methods (Access = protected, Hidden = true) 
        
        
        function [regression_type] = get_regression_type(self)
            regression_type = self.user_inputs.tab_2_lr.regression_type;
            if ~ismember(regression_type, [1 2 3 4 5 6])
                error(['Value error for regression type. Should be integer between 1 and 6.']);
            end
        end
        
        
        function [iterations] = get_iterations(self)
            iterations = self.user_inputs.tab_2_lr.iterations;
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
        
        
        function [burnin] = get_burnin(self)
            burnin = self.user_inputs.tab_2_lr.burnin;
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
        
        
        function [model_credibility] =  get_model_credibility(self)
            model_credibility = self.user_inputs.tab_2_lr.model_credibility;
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
        
        
        function [b] = get_b(self)
            b = self.user_inputs.tab_2_lr.b;
            if ~ismember(class(b), ["char" "double"])
                error(['Type error for b. Should be float or scalar array.']);
            end
            if ischar(b)
                b = iu.char_to_array(b);
                if any(isnan(str2double(b)))
                    error(['Type error for b. All elements should be scalars.']);
                else
                    b = str2double(b)';
                end
            end
            if numel(b) ~= numel(self.exogenous_variables) && numel(b) ~= 1
                error(['Dimension error for b. Dimension of b and exogenous don''t match.']);
            end
            if any(isnan(b)) || any(isinf(b))
                error(['Type error for b. All elements should be scalars.']);
            end
        end
        
        
        function [V] = get_V(self)
            V = self.user_inputs.tab_2_lr.V;
            if ~ismember(class(V), ["char" "double"])
                error(['Type error for V. Should be float or scalar array.']);
            end
            if ischar(V)
                V = iu.char_to_array(V);
                if any(isnan(str2double(V)))
                    error(['Type error for V. All elements should be scalars.']);
                else
                    V = str2double(V)';
                end
            end
            if numel(V) ~= numel(self.exogenous_variables) && numel(V) ~= 1
                error(['Dimension error for V. Dimension of V and exogenous don''t match.']);
            end
            if any(isnan(V)) || any(isinf(V))
                error(['Type error for V. All elements should be scalars.']);
            end
            if any(V < 0)
                error(['Value error for V. Should be positive scalars.']);
            end
        end        
        
        
        function [alpha] = get_alpha(self)
            alpha = self.user_inputs.tab_2_lr.alpha;
            if ~ismember(class(alpha), ["char" "double"])
                error(['Type error for alpha. Should be float or integer.']);
            end
            if ischar(alpha)
                if isnan(str2double(alpha))
                    error(['Type error for alpha. Should be float or integer.']);
                else
                    alpha = str2double(alpha);
                end
            end
            if alpha <= 0
                error(['Value error for alpha. Should be strictly positive.']);
            end
        end         
        
        
        function [delta] = get_delta(self)
            delta = self.user_inputs.tab_2_lr.delta;
            if ~ismember(class(delta), ["char" "double"])
                error(['Type error for delta. Should be float or integer.']);
            end
            if ischar(delta)
                if isnan(str2double(delta))
                    error(['Type error for delta. Should be float or integer.']);
                else
                    delta = str2double(delta);
                end
            end
            if delta <= 0
                error(['Value error for delta. Should be strictly positive.']);
            end
        end        
       
        
        function [g] = get_g(self)
            g = self.user_inputs.tab_2_lr.g;
            if ~ismember(class(g), ["char" "double"])
                error(['Type error for g. Should be float or scalar array.']);
            end
            if ischar(g)
                g = iu.char_to_array(g);
                if any(isnan(str2double(g)))
                    error(['Type error for g. All elements should be scalars.']);
                else
                    g = str2double(g)';
                end
            end
            if any(isnan(g)) || any(isinf(g))
                error(['Type error for g. All elements should be scalars.']);
            end
        end        
        
        
        function [Q] = get_Q(self)
            Q = self.user_inputs.tab_2_lr.Q;
            if ~ismember(class(Q), ["char" "double"])
                error(['Type error for Q. Should be float or scalar array.']);
            end
            if ischar(Q)
                Q = iu.char_to_array(Q);
                if any(isnan(str2double(Q)))
                    error(['Type error for Q. All elements should be scalars.']);
                else
                    Q = str2double(Q)';
                end
            end
            if any(isnan(Q)) || any(isinf(Q))
                error(['Type error for Q. All elements should be scalars.']);
            end
            if any(Q < 0)
                error(['Value error for Q. Should be positive scalars.']);
            end
        end          
        
        
        function [tau] = get_tau(self)
            tau = self.user_inputs.tab_2_lr.tau;
            if ~ismember(class(tau), ["char" "double"])
                error(['Type error for tau. Should be float or integer.']);
            end
            if ischar(tau)
                if isnan(str2double(tau))
                    error(['Type error for tau. Should be float or integer.']);
                else
                    tau = str2double(tau);
                end
            end
            if tau <= 0
                error(['Value error for tau. Should be strictly positive.']);
            end
        end             
        
        
        function [thinning] = get_thinning(self)
            thinning = self.user_inputs.tab_2_lr.thinning;
            if ~islogical(thinning)
                error(['Type error for thinning. Should be boolean.']);
            end
        end        

        
        function [thinning_frequency] = get_thinning_frequency(self)
            thinning_frequency = self.user_inputs.tab_2_lr.thinning_frequency;
            if ~(ischar(thinning_frequency) || iu.is_integer(thinning_frequency))
                error(['Type error for thinning frequency. Should be integer.']);
            end
            if ~iu.is_empty(thinning_frequency) && ischar(thinning_frequency)
                if iu.is_digit(thinning_frequency)
                    thinning_frequency = str2double(thinning_frequency);
                else
                    error(['Type error for thinning frequency. Should be positive integer.']);
                end
            end  
            if iu.is_integer(thinning_frequency) && thinning_frequency <= 0
                error(['Value error for thinning frequency. Should be positive integer.']);
            end            
        end         
        
        
        function [Z_variables] = get_Z_variables(self)
            Z_variables = self.user_inputs.tab_2_lr.Z_variables;
            if self.regression_type == 5
                if iu.is_empty(Z_variables) || ~ismember(class(Z_variables), ["string" "char"])
                    error(['Type error for Z variables. Should be string or char.']);
                end            
                Z_variables = iu.char_to_array(Z_variables);
            else
                Z_variables = " ";
            end
        end        

        
        function [q] = get_q(self)
            q = self.user_inputs.tab_2_lr.q;
            if ~(ischar(q) || iu.is_integer(q))
                error(['Type error for q. Should be integer.']);
            end
            if ~iu.is_empty(q) && ischar(q)
                if iu.is_digit(q)
                    q = str2double(q);
                else
                    error(['Type error for q. Should be positive integer.']);
                end
            end  
            if iu.is_integer(q) && q <= 0
                error(['Value error for q. Should be positive integer.']);
            end            
        end          
        
        
        function [p] = get_p(self)
            p = self.user_inputs.tab_2_lr.p;
            if ~ismember(class(p), ["char" "double"])
                error(['Type error for p. Should be float or scalar array.']);
            end
            if ischar(p)
                p = iu.char_to_array(p);
                if any(isnan(str2double(p)))
                    error(['Type error for p. All elements should be scalars.']);
                else
                    p = str2double(p)';
                end
            end
            if size(p,1) ~= self.q && size(p,1) ~= 1
                error(['Dimension error for p. Dimension of p and lag length q don''t match.']);
            end
            if any(isnan(p)) || any(isinf(p))
                error(['Type error for p. All elements should be scalars.']);
            end
        end         
        

        function [H] = get_H(self)
            H = self.user_inputs.tab_2_lr.H;
            if ~ismember(class(H), ["char" "double"])
                error(['Type error for H. Should be float or scalar array.']);
            end
            if ischar(H)
                H = iu.char_to_array(H);
                if any(isnan(str2double(H)))
                    error(['Type error for H. All elements should be scalars.']);
                else
                    H = str2double(H)';
                end
            end
            if size(H,1) ~= self.q && size(H,1) ~= 1
                error(['Dimension error for H. Dimension of H and lag length q don''t match.']);
            end            
            if any(isnan(H)) || any(isinf(H))
                error(['Type error for H. All elements should be scalars.']);
            end
            if any(H < 0)
                error(['Value error for H. Should be positive scalars.']);
            end
        end
        

        function [constant] = get_constant(self)
            constant = self.user_inputs.tab_2_lr.constant;
            if ~islogical(constant)
                error(['Type error for constant. Should be boolean.']);
            end
        end 
        
        
        function [b_constant] = get_b_constant(self)
            b_constant = self.user_inputs.tab_2_lr.b_constant;
            if ~ismember(class(b_constant), ["char" "double"])
                error(['Type error for b (constant). Should be float or integer.']);
            end
            if ischar(b_constant)
                if isnan(str2double(b_constant))
                    error(['Type error for b (constant). Should be float or integer.']);
                else
                    b_constant = str2double(b_constant);
                end
            end
        end   
        
        
        function [V_constant] = get_V_constant(self)
            V_constant = self.user_inputs.tab_2_lr.V_constant;
            if ~ismember(class(V_constant), ["char" "double"])
                error(['Type error for V (constant). Should be float or integer.']);
            end
            if ischar(V_constant)
                if isnan(str2double(V_constant))
                    error(['Type error for V (constant). Should be float or integer.']);
                else
                    V_constant = str2double(V_constant);
                end
            end
            if V_constant <= 0
                error(['Value error for V (constant). Should be strictly positive.']);
            end
        end         
        
        
        function [trend] = get_trend(self)
            trend = self.user_inputs.tab_2_lr.trend;
            if ~islogical(trend)
                error(['Type error for trend. Should be boolean.']);
            end
        end 
        
        
        function [b_trend] = get_b_trend(self)
            b_trend = self.user_inputs.tab_2_lr.b_trend;
            if ~ismember(class(b_trend), ["char" "double"])
                error(['Type error for b (trend). Should be float or integer.']);
            end
            if ischar(b_trend)
                if isnan(str2double(b_trend))
                    error(['Type error for b (trend). Should be float or integer.']);
                else
                    b_trend = str2double(b_trend);
                end
            end
        end   
        
        
        function [V_trend] = get_V_trend(self)
            V_trend = self.user_inputs.tab_2_lr.V_trend;
            if ~ismember(class(V_trend), ["char" "double"])
                error(['Type error for V (trend). Should be float or integer.']);
            end
            if ischar(V_trend)
                if isnan(str2double(V_trend))
                    error(['Type error for V (trend). Should be float or integer.']);
                else
                    V_trend = str2double(V_trend);
                end
            end
            if V_trend <= 0
                error(['Value error for V (trend). Should be strictly positive.']);
            end
        end           
        
        
        function [quadratic_trend] = get_quadratic_trend(self)
            quadratic_trend = self.user_inputs.tab_2_lr.quadratic_trend;
            if ~islogical(quadratic_trend)
                error(['Type error for quadratic trend. Should be boolean.']);
            end
        end 
        
        
        function [b_quadratic_trend] = get_b_quadratic_trend(self)
            b_quadratic_trend = self.user_inputs.tab_2_lr.b_quadratic_trend;
            if ~ismember(class(b_quadratic_trend), ["char" "double"])
                error(['Type error for b (quadratic trend). Should be float or integer.']);
            end
            if ischar(b_quadratic_trend)
                if isnan(str2double(b_quadratic_trend))
                    error(['Type error for b (quadratic trend). Should be float or integer.']);
                else
                    b_quadratic_trend = str2double(b_quadratic_trend);
                end
            end
        end   
        
        
        function [V_quadratic_trend] = get_V_quadratic_trend(self)
            V_quadratic_trend = self.user_inputs.tab_2_lr.V_quadratic_trend;
            if ~ismember(class(V_quadratic_trend), ["char" "double"])
                error(['Type error for V (quadratic trend). Should be float or integer.']);
            end
            if ischar(V_quadratic_trend)
                if isnan(str2double(V_quadratic_trend))
                    error(['Type error for V (quadratic trend). Should be float or integer.']);
                else
                    V_quadratic_trend = str2double(V_quadratic_trend);
                end
            end
            if V_quadratic_trend <= 0
                error(['Value error for V (quadratic trend). Should be strictly positive.']);
            end
        end         
        
        
        function [insample_fit] = get_insample_fit(self)
            insample_fit = self.user_inputs.tab_2_lr.insample_fit;
            if ~islogical(insample_fit)
                error(['Type error for in-sample fit. Should be boolean.']);
            end
        end         
        
        
        function [marginal_likelihood] = get_marginal_likelihood(self)
            marginal_likelihood = self.user_inputs.tab_2_lr.marginal_likelihood;
            if ~islogical(marginal_likelihood)
                error(['Type error for marginal likelihood. Should be boolean.']);
            end
        end           
        
        
        function [hyperparameter_optimization] = get_hyperparameter_optimization(self)
            hyperparameter_optimization = self.user_inputs.tab_2_lr.hyperparameter_optimization;
            if ~islogical(hyperparameter_optimization)
                error(['Type error for hyperparameter optimization. Should be boolean.']);
            end
        end         
        
        
        function [optimization_type] = get_optimization_type(self)
            optimization_type = self.user_inputs.tab_2_lr.optimization_type;
            if ~ismember(optimization_type, [1 2])
                error(['Value error for optimization type. Should be 1 or 2.']);
            end
        end
        
        
        function [endogenous, exogenous, dates] = get_insample_data(self)
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
        
        
        function [Z] = get_heteroscedastic_data(self)
            % load data only if specified model is heteroscedastic regression
            if self.regression_type == 5
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.data_file);
                % then load data file
                data = iu.load_data(self.project_path, self.data_file);
                % check that Z variables are found in data
                iu.check_variables(data, self.data_file, self.Z_variables, 'Z variables');
                % check that the start and end dates can be found in the file        
                iu.check_dates(data, self.data_file, self.start_date, self.end_date);
                % recover Z variables
                [Z] = iu.fetch_data(data, self.data_file, self.start_date, ...
                         self.end_date, self.Z_variables, 'Z variables');   
                % if dimensions are not consistent with g or Q, raise error
                h = size(Z,2);
                if (numel(self.g)~=1 && h ~= numel(self.g)) || (numel(self.Q)~=1 && h ~= numel(self.Q))
                    error(['Dimension error for heteroscedastic regressors Z. The dimensions of g, Q and Z are not consistent.']);
                end
            % if model is not heteroscedastic regression, return empty list
            else
                Z = [];
            end
        end        


        function [X_p, y_p, Z_p, forecast_dates] = get_forecast_data(self)
            % if forecast is selected, load data
            if self.forecast
                % check that data path and files are valid
                iu.check_file_path(self.project_path, self.forecast_file);
                % then load data file
                data = iu.load_data(self.project_path, self.forecast_file);
                % define the number of forecast periods
                periods = size(data,1);
                % check that exogenous variables are found in data
                iu.check_variables(data, self.forecast_file, self.exogenous_variables, 'Exogenous variables');
                % recover endogenous and exogenous data
                y_p = iu.fetch_forecast_data(data, [], self.endogenous_variables, ...
                self.forecast_file, self.forecast_evaluation, periods, 'Endogenous variable'); 
                X_p = iu.fetch_forecast_data(data, [], self.exogenous_variables, ...
                self.forecast_file, true, periods, 'Exogenous variable');
                % if model is heteroscedastic, Z variables must be obtained as well
                if self.regression_type == 5
                    iu.check_variables(data, self.forecast_file, self.Z_variables, 'Z variables');
                    % recover Z variables
                    Z_p = iu.fetch_forecast_data(data, [], self.Z_variables, ...
                    self.forecast_file, true, periods, 'Z variable');
                else
                    Z_p = [];
                end
                % recover forecast dates
                end_date = self.dates(end);
                forecast_dates = iu.generate_forecast_dates(end_date, periods, self.frequency);
            % if forecasts is not selected, return empty data
            else
                X_p = [];
                y_p = [];
                Z_p = [];
                forecast_dates = [];                
            end
        end
        
    end
    
end
        
        