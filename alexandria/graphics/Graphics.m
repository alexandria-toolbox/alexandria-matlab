classdef Graphics < handle & RegressionGraphics & VectorAutoregressionGraphics
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        model
        path
        clear_folder
    end    

    properties (GetAccess = public, SetAccess = {?RegressionGraphics, ?VectorAutoregressionGraphics})
        complementary_information
    end


    %---------------------------------------------------
    % Methods
    %---------------------------------------------------

    
    methods (Access = public) 

   
        function self = Graphics(model, varargin)

            % constructor for the Graphics class
            
            % allow for optional arguments
            parser = inputParser;
            default_complementary_information = struct;
            default_path = [];
            default_clear_folder = false;
            addRequired(parser, 'model');
            addParameter(parser, 'complementary_information', default_complementary_information);
            addParameter(parser, 'path', default_path);
            addParameter(parser, 'clear_folder', default_clear_folder);
            parse(parser, model, varargin{:});
            self.model = model;
            self.complementary_information = parser.Results.complementary_information;
            self.path = parser.Results.path;
            self.clear_folder = parser.Results.clear_folder;
            % initialize save folder
            self.initialize_folder();
            % complement information with possible missing elements
            self.complete_information();
        end  
        

        function insample_fit_graphics(self, show, save)
            model_class = self.complementary_information.model_class;
            % if model is linear regression, make regression insample graphics
            if model_class == 1
                self.regression_fitted(show, save);
                self.regression_residuals(show, save);
            % if model is vector autoregression, make VAR insample graphics
            elseif model_class == 2 || model_class == 3
                self.var_fitted(show, save);
                self.var_residuals(show, save);
                self.var_shocks(show, save);
                self.var_steady_state(show, save);
            end
        end


        function forecast_graphics(self, show, save)
            model_class = self.complementary_information.model_class;
            % if model is linear regression, make regression forecast graphics
            if model_class == 1
                self.regression_forecasts(show, save);
            % if model is vector autoregression, make VAR forecast graphics
            elseif model_class == 2 || model_class == 3
                self.var_forecasts(show, save);
            end
        end


        function conditional_forecast_graphics(self, show, save)
            model_class = self.complementary_information.model_class;
            % if model is vector autoregression, make VAR forecast graphics
            if model_class == 2 || model_class == 3
                self.var_conditional_forecasts(show, save);
            end
        end


        function irf_graphics(self, show, save)
            model_class = self.complementary_information.model_class;
            % if model is vector autoregression, make VAR IRF graphics
            if model_class == 2 || model_class == 3
                self.var_irf(show, save);
            end
        end


        function fevd_graphics(self, show, save)
            model_class = self.complementary_information.model_class;
            % if model is vector autoregression, make VAR FEVD graphics
            if model_class == 2 || model_class == 3
                self.var_fevd(show, save);
            end
        end


        function hd_graphics(self, show, save)
            model_class = self.complementary_information.model_class;
            % if model is vector autoregression, make VAR HD graphics
            if model_class == 2 || model_class == 3
                self.var_hd(show, save);
            end
        end

    end
    

    methods (Access = protected, Hidden = true)
        

        function  initialize_folder(self)

            % check if path to save folder is defined, otherwise default is current directory
            if isempty(self.path)
                self.path = pwd;
            end
            % clear save folder and re-initialize if activated
            if self.clear_folder && exist(self.path) == 7
                rmdir(self.path, 's');
            end
            % create path if it does not exist
            if ~(exist(self.path) == 7)
                mkdir(self.path);
            end
        end


        function complete_information(self)
            % add general model information
            self.complete_model_information();      
            % if model is linear regression, add regression elements
            if self.complementary_information.model_class == 1
                self.complete_regression_information();        
            % if model is vecto autoregression, add var elements
            elseif self.complementary_information.model_class == 2
                self.complete_var_information();
            % if model is VEC/VARMA, add var elements (vec and varma just recycle VAR functions)
            elseif self.complementary_information.model_class == 3
                self.complete_var_information();               
            end
        end


        function complete_model_information(self)
            % recover and add common model elements
            [model_name model_class model_type] = iu.identify_model(self.model);
            self.complementary_information.model_name = model_name;
            self.complementary_information.model_class = model_class;
            self.complementary_information.model_type = model_type;
        end


    end
    
end
