classdef InputProcessor < handle & RegressionProcessor & VectorAutoregressionProcessor ...
                          & VecVarmaProcessor
    
    
    properties (GetAccess = public, SetAccess= protected)
        user_inputs
        model
        endogenous_variables
        exogenous_variables
        frequency
        start_date
        end_date
        project_path
        data_file
        progress_bar
        create_graphics
        save_results
        forecast
        forecast_credibility
        conditional_forecast
        conditional_forecast_credibility
        irf
        irf_credibility
        fevd
        fevd_credibility
        hd
        hd_credibility
        forecast_periods
        conditional_forecast_type
        forecast_file
        conditional_forecast_file
        forecast_evaluation
        irf_periods
        structural_identification
        structural_identification_file
        estimation_start
        estimation_end
    end    


    properties (GetAccess = public, SetAccess = {?RegressionProcessor, ...
                ?VectorAutoregressionProcessor, ?VecVarmaProcessor})
        results_information
        graphics_information
    end

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------    
    
    
    methods (Access = public)    
    
    
        function self = InputProcessor(user_inputs)    
            % save inputs as attributes
            self.user_inputs = user_inputs;         
        end
        

        function process_input(self)   
            % get inputs of tab 1 of graphical interface (common to all models)
            self.tab_1_inputs();
            % get inputs from tab 2 of graphical interface (model-specific)
            self.tab_2_inputs();
            % then get inputs of tab 3 of graphical interface (again common to all models)
            self.tab_3_inputs();
            % finally, load all relevant data files
            self.data_inputs();
        end


        function input_timer(self, tag)
            if isequal(tag, 'start')
                self.estimation_start = now;
            elseif isequal(tag, 'end')
                self.estimation_end = now;
            end
        end


        function make_results_information(self)
            % initialize information dictionary
            self.results_information = struct;
            % make general information
            self.make_general_information();
            % if model is model 1, additionally make information for linear regression
            if self.model == 1
                self.make_regression_information();
            % if model is model 2, additionally make information for VAR models
            elseif self.model == 2
                self.make_var_information();
            % if model is model 3, additionally make information for VEC/VARMA models
            elseif self.model == 3
                self.make_vec_varma_information();                
            end
            % finally add complementary information for applications
            self.make_application_information();
        end


        function make_graphics_information(self)
            % initialize graphics dictionary
            self.graphics_information = struct;
            % get endogenous variable
            self.graphics_information.endogenous_variables = self.endogenous_variables;
            % get exogenous variables
            if isempty(self.exogenous_variables)
                self.graphics_information.exogenous_variables = [];
            else
                self.graphics_information.exogenous_variables = self.exogenous_variables;
            end
            % if model is model 1, additionally make information for linear regression
            if self.model == 1
                self.make_regression_graphics_information();
            % if model is model 2, additionally make information for vector autoregression
            elseif self.model == 2
                self.make_var_graphics_information();
            % if model is model 3, additionally make information for VEC/VARMA
            elseif self.model == 3
                self.make_vec_varma_graphics_information();                
            end
        end


    end

    
    %---------------------------------------------------
    % Methods (Access = private)
    %---------------------------------------------------    
    
    
    methods (Access = protected, Hidden = true)    
    
    
        function tab_1_inputs(self)
            % recover model
            self.model = self.get_model();
            % get list of endogenous variables
            self.endogenous_variables = self.get_endogenous_variables();
            % get list of exogenous variables
            self.exogenous_variables = self.get_exogenous_variables();
            % get data frequency
            self.frequency = self.get_frequency();
            % get sample periods
            [self.start_date, self.end_date] = self.get_sample_dates();
            % get path to project folder
            self.project_path = self.get_project_path();
            % get name of data file
            self.data_file = self.get_data_file();
            % get progress bar decision
            self.progress_bar = self.get_progress_bar();
            % get graphic creation decision
            self.create_graphics = self.get_create_graphics();
            % get result saving decision
            self.save_results = self.get_save_results();
        end
        
        
        function tab_2_inputs(self)
            % if model is model 1, get user inputs for linear regression
            if self.model == 1
                self.regression_inputs();
            % if model is model 2, get user inputs for vector autoregression
            elseif self.model == 2
                self.vector_autoregression_inputs();
            % if model is model 3, get user inputs for vec/varma
            elseif self.model == 3
                self.vec_varma_inputs();                
            end
        end
            

        function tab_3_inputs(self)
            % recover forecast decision
            self.forecast = self.get_forecast();
            % recover forecast credibility interval
            self.forecast_credibility = self.get_forecast_credibility();
            % recover conditional forecast decision
            self.conditional_forecast = self.get_conditional_forecast();
            % recover conditional forecast credibility interval
            self.conditional_forecast_credibility = self.get_conditional_forecast_credibility();      
            % recover irf decision
            self.irf = self.get_irf();
            % recover irf credibility interval
            self.irf_credibility = self.get_irf_credibility();
            % recover fevd decision
            self.fevd = self.get_fevd();
            % recover fevd credibility interval
            self.fevd_credibility = self.get_fevd_credibility();
            % recover hd decision
            self.hd = self.get_hd();
            % recover hd credibility interval
            self.hd_credibility = self.get_hd_credibility();
            % recover number of forecast periods
            self.forecast_periods = self.get_forecast_periods();
            % recover type of conditional forecast
            self.conditional_forecast_type = self.get_conditional_forecast_type();
            % recover name of forecast input file
            self.forecast_file = self.get_forecast_file();
            % recover name of conditional forecast input file
            self.conditional_forecast_file = self.get_conditional_forecast_file();
            % recover forecast evaluation decision
            self.forecast_evaluation = self.get_forecast_evaluation();
            % recover number of irf periods
            self.irf_periods = self.get_irf_periods();     
            % recover type of structural identification
            self.structural_identification = self.get_structural_identification();
            % recover type of structural identification
            self.structural_identification_file = self.get_structural_identification_file();
        end
        
        
        function data_inputs(self)
            % if model is model 1, get data for linear regression 
            if self.model == 1
                self.regression_data();
            % else, if model is model 2, get data for vector autoregression
            elseif self.model == 2
                self.vector_autoregression_data();
            % else, if model is model 3, get data for vec/varma
            elseif self.model == 3
                self.vec_varma_data();
            end
        end
        
        
        function [model] = get_model(self)
        model = self.user_inputs.tab_1.model;
            if ~ismember(model, [1 2 3])
                error(['Value error for model. Should be 1, 2 or 3.']);
            end
        end
        

        function [endogenous_variables] = get_endogenous_variables(self)
            endogenous_variables = self.user_inputs.tab_1.endogenous_variables;
            if iu.is_empty(endogenous_variables) || ...
               ~ismember(class(endogenous_variables), ["string" "char"])
                error(['Type error for endogenous variables. Should non-empty char or string.']);
            end
            endogenous_variables = iu.char_to_array(endogenous_variables);
            if self.model == 1 && numel(endogenous_variables) ~= 1
                error(['Dimension error for endogenous variable. Linear regression model shoud specify one and only one endogenous variable.']);
            end
        end
            
            
        function [exogenous_variables] = get_exogenous_variables(self)
            exogenous_variables = self.user_inputs.tab_1.exogenous_variables;
            if ~ismember(class(exogenous_variables), ["string" "char"])
                error(['Type error for exogenous variables. Should non-empty char or string.']);
            end
            exogenous_variables = iu.char_to_array(exogenous_variables);
        end
        
        
        function [frequency] = get_frequency(self)
            frequency = self.user_inputs.tab_1.frequency;
            if ~ismember(frequency, [1 2 3 4 5 6])
                error(['Value error for frequency. Should be int between 1 and 6.']);
            end
        end
        
        
        function [start_date, end_date] = get_sample_dates(self)
            sample_dates = self.user_inputs.tab_1.sample;
            if iu.is_empty(sample_dates) || ~ismember(class(sample_dates), ["string" "char"])
                error(['Type error for sample dates. Should non-empty list of strings.']);
            end
            sample_dates = iu.char_to_array(sample_dates);
            start_date = char(sample_dates(1));
            end_date = char(sample_dates(2));
        end
        
        
        function [project_path] = get_project_path(self)
            project_path = self.user_inputs.tab_1.project_path;
            if ~ismember(class(project_path), ["string" "char"])
                error(['Type error for project folder path. Should be char.']);
            end
            project_path = iu.fix_char(project_path);
            if iu.is_empty(project_path)
                error(['Value error for project folder path. Should be non-empty char.']);
            end
        end
        
        
        function [data_file] = get_data_file(self)
            data_file = self.user_inputs.tab_1.data_file;
            if ~ismember(class(data_file), ["string" "char"])
                error(['Type error for project data file. Should be char.']);
            end            
            data_file = iu.fix_char(data_file);
            if iu.is_empty(data_file)
                error(['Value error for project data file. Should be non-empty char.']);
            end
        end
        
        
        function [progress_bar] = get_progress_bar(self)
            progress_bar = self.user_inputs.tab_1.progress_bar;
            if ~islogical(progress_bar)
                error(['Type error for progress bar. Should be boolean.']);
            end
        end
        
        
        function [create_graphics] = get_create_graphics(self)
            create_graphics = self.user_inputs.tab_1.create_graphics;
            if ~islogical(create_graphics)
                error(['Type error for create graphics. Should be boolean.']);
            end           
        end
        
        
        function [save_results] = get_save_results(self)
            save_results = self.user_inputs.tab_1.save_results;
            if ~islogical(save_results)
                error(['Type error for save results. Should be boolean.']);
            end           
        end
        
        
        function [forecast] = get_forecast(self)
            forecast = self.user_inputs.tab_3.forecast;
            if ~islogical(forecast)
                error(['Type error for forecasts. Should be boolean.']);
            end
        end
        
        
        function [forecast_credibility] = get_forecast_credibility(self)
            forecast_credibility = self.user_inputs.tab_3.forecast_credibility;
            if ~ismember(class(forecast_credibility), ["char" "double"])
                error(['Type error for forecasts credibility level. Should be float between 0 and 1.']);
            end
            if ischar(forecast_credibility)
                if isnan(str2double(forecast_credibility))
                    error(['Type error for forecasts credibility level. Should be float between 0 and 1.']);
                else
                    forecast_credibility = str2double(forecast_credibility);
                end
            end
            if forecast_credibility <= 0 || forecast_credibility >= 1
                error(['Value error for forecasts credibility level. Should be float between 0 and 1 (not included).']);
            end
        end
                    
        
        function [conditional_forecast] = get_conditional_forecast(self)
            conditional_forecast = self.user_inputs.tab_3.conditional_forecast;
            if ~islogical(conditional_forecast)
                error(['Type error for conditional forecasts. Should be boolean.']);
            end
        end        
        

        function [conditional_forecast_credibility] = get_conditional_forecast_credibility(self)
            conditional_forecast_credibility = self.user_inputs.tab_3.conditional_forecast_credibility;
            if ~ismember(class(conditional_forecast_credibility), ["char" "double"])
                error(['Type error for conditional forecasts credibility level. Should be float between 0 and 1.']);
            end
            if ischar(conditional_forecast_credibility)
                if isnan(str2double(conditional_forecast_credibility))
                    error(['Type error for conditional forecasts credibility level. Should be float between 0 and 1.']);
                else
                    conditional_forecast_credibility = str2double(conditional_forecast_credibility);
                end
            end
            if conditional_forecast_credibility <= 0 || conditional_forecast_credibility >= 1
                error(['Value error for conditional forecasts credibility level. Should be float between 0 and 1 (not included).']);
            end
        end
        
        
        function [irf] = get_irf(self)
            irf = self.user_inputs.tab_3.irf;
            if ~islogical(irf)
                error(['Type error for impulse response function. Should be boolean.']);
            end
        end            
        
        
        function [irf_credibility] = get_irf_credibility(self)
            irf_credibility = self.user_inputs.tab_3.irf_credibility;
            if ~ismember(class(irf_credibility), ["char" "double"])
                error(['Type error for irf credibility level. Should be float between 0 and 1.']);
            end
            if ischar(irf_credibility)
                if isnan(str2double(irf_credibility))
                    error(['Type error for irf credibility level. Should be float between 0 and 1.']);
                else
                    irf_credibility = str2double(irf_credibility);
                end
            end
            if irf_credibility <= 0 || irf_credibility >= 1
                error(['Value error for irf credibility level. Should be float between 0 and 1 (not included).']);
            end
        end    
        
        
        function [fevd] = get_fevd(self)
            fevd = self.user_inputs.tab_3.fevd;
            if ~islogical(fevd)
                error(['Type error for forecast error variance decomposition. Should be boolean.']);
            end
        end   
        
        
        function [fevd_credibility] =  get_fevd_credibility(self)
            fevd_credibility = self.user_inputs.tab_3.fevd_credibility;
            if ~ismember(class(fevd_credibility), ["char" "double"])
                error(['Type error for fevd credibility level. Should be float between 0 and 1.']);
            end
            if ischar(fevd_credibility)
                if isnan(str2double(fevd_credibility))
                    error(['Type error for fevd credibility level. Should be float between 0 and 1.']);
                else
                    fevd_credibility = str2double(fevd_credibility);
                end
            end
            if fevd_credibility <= 0 || fevd_credibility >= 1
                error(['Value error for fevd credibility level. Should be float between 0 and 1 (not included).']);
            end
        end 
        
        
        function [hd] = get_hd(self)
            hd = self.user_inputs.tab_3.hd;
            if ~islogical(hd)
                error(['Type error for historical decomposition. Should be boolean.']);
            end
        end 
        
        
        function [hd_credibility] =  get_hd_credibility(self)
            hd_credibility = self.user_inputs.tab_3.hd_credibility;
            if ~ismember(class(hd_credibility), ["char" "double"])
                error(['Type error for historical decomposition credibility level. Should be float between 0 and 1.']);
            end
            if ischar(hd_credibility)
                if isnan(str2double(hd_credibility))
                    error(['Type error for historical decomposition credibility level. Should be float between 0 and 1.']);
                else
                    hd_credibility = str2double(hd_credibility);
                end
            end
            if hd_credibility <= 0 || hd_credibility >= 1
                error(['Value error for historical decomposition credibility level. Should be float between 0 and 1 (not included).']);
            end
        end 
        
        
        function [forecast_periods] = get_forecast_periods(self)
            forecast_periods = self.user_inputs.tab_3.forecast_periods;
            if ~(ischar(forecast_periods) || iu.is_integer(forecast_periods))
                error(['Type error for forecast periods. Should be integer.']);
            end
            if ismember(self.model, [2]) && (self.forecast || self.conditional_forecast) && isempty(forecast_periods)
                error(['Type error for forecast periods. Should be integer.']);
            end
            if ~iu.is_empty(forecast_periods) && ischar(forecast_periods)
                if iu.is_digit(forecast_periods)
                    forecast_periods = str2double(forecast_periods);
                else
                    error(['Type error for forecast periods. Should be positive integer.']);
                end
            end
            if iu.is_integer(forecast_periods) && forecast_periods <= 0
                error(['Value error for forecast periods. Should be positive integer.']);
            end
        end
        
        
        function [conditional_forecast_type] = get_conditional_forecast_type(self)
            conditional_forecast_type = self.user_inputs.tab_3.conditional_forecast_type;
            if ~ismember(conditional_forecast_type, [1, 2])
                error(['Value error for conditional forecast type. Should be 1 or 2.']);
            end
        end
        
        
        function [forecast_file] = get_forecast_file(self)
            forecast_file = self.user_inputs.tab_3.forecast_file;
            if ~ischar(forecast_file)
                error(['Type error for forecast file. Should be char.']);
            end
            forecast_file = iu.fix_char(forecast_file);
        end
        
        
        function [conditional_forecast_file] = get_conditional_forecast_file(self)
            conditional_forecast_file = self.user_inputs.tab_3.conditional_forecast_file;
            if ~ischar(conditional_forecast_file)
                error(['Type error for conditional forecast file. Should be char.']);
            end
            conditional_forecast_file = iu.fix_char(conditional_forecast_file);
        end
        
        
        function [forecast_evaluation] = get_forecast_evaluation(self)
            forecast_evaluation = self.user_inputs.tab_3.forecast_evaluation;
            if ~islogical(forecast_evaluation)
                error(['Type error for forecast evaluation. Should be boolean.']);
            end
        end
        
        
        function [irf_periods] = get_irf_periods(self)
            irf_periods = self.user_inputs.tab_3.irf_periods;
            if ~(ischar(irf_periods) || iu.is_integer(irf_periods))
                error(['Type error for irf periods. Should be integer.']);
            end
            if ismember(self.model, [2]) && self.irf && isempty(irf_periods)
                error(['Type error for irf periods. Should be integer.']);
            end
            if ~iu.is_empty(irf_periods) && ischar(irf_periods)
                if iu.is_digit(irf_periods)
                    irf_periods = str2double(irf_periods);
                else
                    error(['Type error for irf periods. Should be positive integer.']);
                end
            end    
            if iu.is_integer(irf_periods) && irf_periods <= 0
                error(['Value error for irf periods. Should be positive integer.']);
            end
        end
        
        
        function [structural_identification] = get_structural_identification(self)
            structural_identification = self.user_inputs.tab_3.structural_identification;
            if ~ismember(structural_identification, [1 2 3 4])
                error(['Value error for structural identification. Should be 1, 2, 3 or 4.']);
            end
            if self.model == 2 && self.user_inputs.tab_2_var.var_type == 1 && ~ismember(structural_identification, [1 2 3])
                error(['Value error for structural identification. Identification by restriction is not available for maximum likelihood VAR.']);            
            end
            if self.model == 2 && self.user_inputs.tab_2_var.var_type == 7 && ~ismember(structural_identification, [1 4])
                error(['Value error for structural identification. Should be 1 (none) or 4 (restrictions) when selecting a proxy SVAR.']);            
            end
        end
        
        
        function [structural_identification_file] = get_structural_identification_file(self)
            structural_identification_file = self.user_inputs.tab_3.structural_identification_file;
            if ~ischar(structural_identification_file)
                error(['Type error for structural identification file. Should be char.']);
            end
            structural_identification_file = iu.fix_char(structural_identification_file);
        end   


        function make_general_information(self)
            % get estimation start date
            estimation_start = datestr(self.estimation_start, 'yyyy-mm-dd hh:MM:ss');
            self.results_information.estimation_start = estimation_start;
            % get estimation start date
            estimation_end = datestr(self.estimation_end, 'yyyy-mm-dd hh:MM:ss');
            self.results_information.estimation_end = estimation_end;
            % get endogenous variable
            self.results_information.endogenous_variables = self.endogenous_variables;
            % get exogenous variables
            if isequal(self.exogenous_variables,"")
                self.results_information.exogenous_variables = "none";
            else
                self.results_information.exogenous_variables = self.exogenous_variables;
            end
            % get sample start date
            self.results_information.sample_start = self.start_date;
            % get sample end date
            self.results_information.sample_end = self.end_date;
            % get data frequency
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
            self.results_information.frequency = frequency;  
            % get path to project folder
            self.results_information.project_path = self.project_path;
            % get data file
            self.results_information.data_file = self.data_file;
            % get progress bar
            self.results_information.progress_bar = self.progress_bar;
            % get graphics and figures
            self.results_information.create_graphics = self.create_graphics;
            % get result save
            self.results_information.save_results = self.save_results;
        end


        function make_application_information(self)
            % forecasts
            self.results_information.forecast = self.forecast;
            self.results_information.forecast_credibility = self.forecast_credibility;
            % conditional forecasts
            self.results_information.conditional_forecast = self.conditional_forecast;
            self.results_information.conditional_forecast_credibility = self.conditional_forecast_credibility;
            % irf
            self.results_information.irf = self.irf;
            self.results_information.irf_credibility = self.irf_credibility;       
            % fevd
            self.results_information.fevd = self.fevd;
            self.results_information.fevd_credibility = self.fevd_credibility; 
            % fevd
            self.results_information.hd = self.hd;
            self.results_information.hd_credibility = self.hd_credibility; 
            % forecast specifications
            self.results_information.forecast_periods = self.forecast_periods;
            self.results_information.conditional_forecast_type = self.conditional_forecast_type;
            self.results_information.forecast_file = self.forecast_file;
            self.results_information.conditional_forecast_file = self.conditional_forecast_file;
            self.results_information.forecast_evaluation = self.forecast_evaluation;
            % irf specifications
            self.results_information.irf_periods = self.irf_periods;
            self.results_information.structural_identification = self.structural_identification;
            self.results_information.structural_identification_file = self.structural_identification_file;
        end

    end
    
end
