classdef Results < handle & RegressionResults & VectorAutoregressionResults ...
                   & VecVarmaResults
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = {?RegressionResults, ?VectorAutoregressionResults, ...
                                                 ?VecVarmaResults})
        model
        complementary_information
    end    


    properties (GetAccess = public, SetAccess = {?RegressionResults, ?VectorAutoregressionResults, ...
                                                 ?VecVarmaResults})
        input_summary
        estimation_summary
        application_summary
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------

    
    methods (Access = public) 

        
        function self = Results(model, varargin)
            % save attributes
            parser = inputParser;
            default_complementary_information = struct;
            addRequired(parser,'model');
            addParameter(parser, 'complementary_information', default_complementary_information);
            parse(parser, model, varargin{:});
            self.model = model;
            self.complementary_information = parser.Results.complementary_information;
            % complement information with possible missing elements
            self.complete_information();
        end 
        
            
        function make_input_summary(self)
            % initiate string list
            self.input_summary = [];
            % add settings header
            self.add_input_header()
            % add tab 1 settings
            self.add_tab_1_inputs()       
            % add tab 2 settings
            self.add_tab_2_inputs()
            % add tab 3 settings
            self.add_tab_3_inputs()
        end


        function show_input_summary(self)
            % display input summary in console
            cu.print_string_array(self.input_summary);
        end


        function save_input_summary(self, path)
            % check if path exists, and create directory if needed
            cu.check_path(path);
            % generate full path
            full_path = fullfile(path, 'input_summary.txt');
            % write txt file on disk
            input_summary = [cu.alexandria_header();self.input_summary];
            cu.write_string_array(input_summary, full_path);
        end


        function make_estimation_summary(self)
            model_class = self.complementary_information.model_class;
            % if model is linear regression, make regression summary
            if model_class == 1
                self.make_regression_summary();       
            % if model is vector autoregression, make VAR summary
            elseif model_class == 2
                self.make_var_summary();
            % if model is VEC/VARMA, make VAR extension summary
            elseif model_class == 3
                self.make_vec_varma_summary();                
            end
        end


        function show_estimation_summary(self)
            % display estimation summary in console
            cu.print_string_array(self.estimation_summary);
        end

    
        function save_estimation_summary(self, path)
            % check if path exists, and create directory if needed
            cu.check_path(path);
            % generate full path
            full_path = fullfile(path, 'estimation_summary.txt');
            % write txt file on disk
            estimation_summary = [cu.alexandria_header();self.estimation_summary];
            cu.write_string_array(estimation_summary, full_path);
        end


        function make_application_summary(self)
            % initiate application_summary
            self.application_summary = struct;
            % identify model to run relevant application summary
            model_class = self.complementary_information.model_class;
            % if model is linear regression, make regression summary
            if model_class == 1
                self.make_regression_application_summary();
            % if model is vector autoregression, make VAR summary
            elseif model_class == 2
                self.make_var_application_summary();
            % if model is VEC/VARMA, make VAR extension summary
            elseif model_class == 3
                self.make_vec_varma_application_summary();                
            end
        end


        function save_application_summary(self, path)
            % check if path exists, and create directory if needed
            cu.check_path(path);
            % identify model to run relevant application summary
            model_class = self.complementary_information.model_class;
            % if model is linear regression, save regression summary
            if model_class == 1
                self.save_regression_application(path);        
            % if model is vector autoregression, save regression summary
            elseif model_class == 2
                self.save_var_application(path);
            % if model is VEC/VARMA, save VAR extension summary
            elseif model_class == 3
                self.save_vec_varma_application(path);                
            end
        end

    end


    methods (Access = protected, Hidden = true)
        

        function complete_information(self)
            % add general model information
            self.complete_model_information();
            % if model is linear regression, add regression elements
            if self.complementary_information.model_class == 1
                self.complete_regression_information();
            % if model is vector autoregression, add VAR elements
            elseif self.complementary_information.model_class == 2
                self.complete_var_information();
            % if model is VEC/VARMA, add extension elements
            elseif self.complementary_information.model_class == 3
                self.complete_vec_varma_information();                
            end
            % add application information
            self.complete_application_information();
        end


        function complete_model_information(self)
            % recover and add common model elements
            [model_name model_class model_type] = iu.identify_model(self.model);
            self.complementary_information.model_name = model_name;
            self.complementary_information.model_class = model_class;
            self.complementary_information.model_type = model_type;
            if ~isfield(self.complementary_information, 'estimation_start')
                self.complementary_information.estimation_start = '—';
            end
            if ~isfield(self.complementary_information, 'estimation_end')
                self.complementary_information.estimation_end = datestr(now, 'yyyy-mm-dd hh:MM:ss');
            end
            if ~isfield(self.complementary_information, 'sample_start')
                self.complementary_information.sample_start = '';
            end            
            if ~isfield(self.complementary_information, 'sample_end')
                self.complementary_information.sample_end = '';
            end
            if ~isfield(self.complementary_information, 'frequency')
                self.complementary_information.frequency = '—';
            end
            if ~isfield(self.complementary_information, 'project_path')
                self.complementary_information.project_path = '—';
            end
            if ~isfield(self.complementary_information, 'data_file')
                self.complementary_information.data_file = '—';
            end
            if ~isfield(self.complementary_information, 'progress_bar')
                self.complementary_information.progress_bar = '—';
            end
            if ~isfield(self.complementary_information, 'create_graphics')
                self.complementary_information.create_graphics = '—';
            end
            if ~isfield(self.complementary_information, 'save_results')
                self.complementary_information.save_results = '—';
            end
        end

    
        function complete_application_information(self)
            if ~isfield(self.complementary_information, 'forecast')
                self.complementary_information.forecast = '—';
            end
            if ~isfield(self.complementary_information, 'forecast_credibility')
                self.complementary_information.forecast_credibility = '—';
            end            
            if ~isfield(self.complementary_information, 'conditional_forecast')
                self.complementary_information.conditional_forecast = '—';
            end
            if ~isfield(self.complementary_information, 'conditional_forecast_credibility')
                self.complementary_information.conditional_forecast_credibility = '—';
            end
            if ~isfield(self.complementary_information, 'irf')
                self.complementary_information.irf = '—';
            end
            if ~isfield(self.complementary_information, 'irf_credibility')
                self.complementary_information.irf_credibility = '—';
            end
            if ~isfield(self.complementary_information, 'fevd')
                self.complementary_information.fevd = '—';
            end
            if ~isfield(self.complementary_information, 'fevd_credibility')
                self.complementary_information.fevd_credibility = '—';
            end
            if ~isfield(self.complementary_information, 'hd')
                self.complementary_information.hd = '—';
            end
            if ~isfield(self.complementary_information, 'hd_credibility')
                self.complementary_information.hd_credibility = '—';
            end
            if ~isfield(self.complementary_information, 'forecast_periods')
                self.complementary_information.forecast_periods = '—';
            end
            if ~isfield(self.complementary_information, 'conditional_forecast_type')
                self.complementary_information.conditional_forecast_type = '—';
            end
            if ~isfield(self.complementary_information, 'forecast_file')
                self.complementary_information.forecast_file = '—';
            end
            if ~isfield(self.complementary_information, 'conditional_forecast_file')
                self.complementary_information.conditional_forecast_file = '—';
            end
            if ~isfield(self.complementary_information, 'forecast_evaluation')
                self.complementary_information.forecast_evaluation = '—';
            end
            if ~isfield(self.complementary_information, 'irf_periods')
                self.complementary_information.irf_periods = '—';
            end
            if ~isfield(self.complementary_information, 'structural_identification')
                self.complementary_information.structural_identification = '—';
            end
            if ~isfield(self.complementary_information, 'structural_identification_file')
                self.complementary_information.structural_identification_file = '—';
            end
        end


        function add_input_header(self)
            % Alexandria header and estimation date
            lines = string([]);
            lines = [lines;['Estimation date:  ' self.complementary_information.estimation_end]];
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end


        function add_tab_1_inputs(self)
            % initiate lines
            lines = string([]);
            % header for tab 1
            lines = [lines;'Model'];
            lines = [lines;'---------'];
            lines = [lines;' '];
            % model class
            model_class = self.complementary_information.model_class;
            if model_class == 1
                model = 'linear regression';
            elseif model_class == 2
                model = 'vector autoregression';
            elseif model_class == 3
                model = 'vec / varma';                
            end
            lines = [lines;'selected model: ' model];
            % endogenous variables
            endogenous_variables = iu.array_to_char(self.complementary_information.endogenous_variables);
            lines = [lines;'endogenous variables: ' endogenous_variables];
            % exogenous variables
            exogenous_variables = iu.array_to_char(self.complementary_information.exogenous_variables);
            lines = [lines;'exogenous variables: ' exogenous_variables];
            % data frequency
            frequency = self.complementary_information.frequency;
            lines = [lines;'data frequency: ' frequency];
            % sample dates
            sample_start = self.complementary_information.sample_start;
            sample_end = self.complementary_information.sample_end;
            sample_dates = [sample_start ' ' sample_end];
            lines = [lines;'estimation sample: ' sample_dates];
            % project path and data file
            project_path = self.complementary_information.project_path;
            data_file = self.complementary_information.data_file;
            lines = [lines;'path to project folder: ' project_path];
            lines = [lines;'data file: ' data_file];
            % progress bar, graphics and result saves
            if islogical(self.complementary_information.progress_bar)
                progress_bar = cu.bool_to_string(self.complementary_information.progress_bar);
            else
                progress_bar = self.complementary_information.progress_bar;
            end
            if islogical(self.complementary_information.create_graphics)
                create_graphics = cu.bool_to_string(self.complementary_information.create_graphics);
            else
                create_graphics = self.complementary_information.create_graphics;
            end
            if islogical(self.complementary_information.save_results)
                save_results = cu.bool_to_string(self.complementary_information.save_results);
            else
                save_results = self.complementary_information.save_results;
            end
            lines = [lines;'progress bar: ' progress_bar];
            lines = [lines;'create graphics: ' create_graphics];
            lines = [lines;'save_results: ' save_results];
            lines = [lines;' '];
            lines = [lines;' '];
            self.input_summary = [self.input_summary;lines];
        end


        function add_tab_2_inputs(self)
            % tab 2 elements for linear regression
            if self.complementary_information.model_class == 1
                self.add_regression_tab_2_inputs();
            elseif self.complementary_information.model_class == 2
                self.add_var_tab_2_inputs();
            elseif self.complementary_information.model_class == 3
                self.add_vec_varma_tab_2_inputs();                
            end
        end

           
        function add_tab_3_inputs(self)
            % tab 3 elements for linear regression
            if self.complementary_information.model_class == 1
                self.add_regression_tab_3_inputs();
            elseif self.complementary_information.model_class == 2
                self.add_var_tab_3_inputs();
            elseif self.complementary_information.model_class == 3
                self.add_vec_varma_tab_3_inputs();                
            end
        end

    end
    
end