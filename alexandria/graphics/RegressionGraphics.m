classdef RegressionGraphics < handle & Graphics
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        insample_fit
        forecast
        project_path
        endogenous
        insample_dates
        y_p
        forecast_dates
        actual
        fitted
        residuals
        estimates_forecasts
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    


        function [self] = RegressionGraphics(ip, lr)
            % determine which applications will produce graphics
            self.graphics_information(ip);
            % gather information from input processor 
            self.input_information(ip);
            % then gather information from regression model
            self.regression_information(lr);     
        end
        
        
        function make_graphics(self)
            % delete existing graphics folder, if any
            self.delete_graphics_folder();
            % create graphics for in-sample fit, if selected
            self.insample_fit_graphics();
            % create graphics for forecasts, if selected
            self.forecasts_graphics();      
        end
        
        
    end
    
    
    methods (Access = protected, Hidden = true)
        
        
        function graphics_information(self, ip)
            % input processor information: in-sample fit
            self.insample_fit = ip.insample_fit;
            % input processor information: forecast decision
            self.forecast = ip.forecast;
        end
        
        
        function input_information(self, ip)
            % input processor information: project folder
            self.project_path = ip.project_path;
            % input processor information: endogenous
            self.endogenous = ip.endogenous_variables;
            % data specific to in-sample fit
            if self.insample_fit
                % input processor information: in-sample dates
                self.insample_dates = ip.dates  ;
            end
            % data specific to forecasts
            if self.forecast
                % input processor information: endogenous, actual values for predictions
                self.y_p = ip.y_p;
                % input processor information: forecast dates
                self.forecast_dates = ip.forecast_dates;
            end
        end
        
        
        function regression_information(self, lr)
            % data specific to in-sample fit
            if self.insample_fit
                % regression information: endogenous values
                self.actual = lr.y;
                % regression information: fitted
                self.fitted = lr.estimates_fit;
                % regression information: residual estimates
                self.residuals = lr.estimates_residuals;
            end
            % data specific to forecasts
            if self.forecast
                % regression information: forecast estimates
                self.estimates_forecasts = lr.estimates_forecasts;
            end
        end
        
        
        function insample_fit_graphics(self)
            if self.insample_fit
                actual = self.actual;
                fitted = self.fitted;
                residuals = self.residuals;
                dates = self.insample_dates;
                name = convertStringsToChars(self.endogenous(1));
                path = self.project_path;
                gu.fit_single_variable(actual, fitted, dates, name, path);
                gu.residual_single_variable(residuals, dates, name, path);
            end
        end
        
        
        function forecasts_graphics(self)
            if self.forecast
                forecasts = self.estimates_forecasts;
                y_p = self.y_p;
                dates = self.forecast_dates;
                name = convertStringsToChars(self.endogenous(1));
                path = self.project_path;
                gu.ols_forecasts_single_variable(forecasts, y_p, dates, name, path);
            end
        end

        
    end
    
    
end