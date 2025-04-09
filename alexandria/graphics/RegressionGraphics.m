classdef RegressionGraphics < handle
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    


        function [self] = RegressionGraphics()    
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
                self.complementary_information.dates = (1:n)';
            end   
            % forecast_dates
            if ~isempty(self.model.forecast_estimates)
                if self.complementary_information.model_type == 6 && ...
                   isfield(self.complementary_information, 'forecast_dates')
                    forecast_dates = self.complementary_information.forecast_dates;
                else
                    forecast_dates = (1:size(self.model.forecast_estimates,1))';
                end
                self.complementary_information.forecast_dates = forecast_dates;
            end
            % actual
            if ~isfield(self.complementary_information, 'y_p')
                self.complementary_information.y_p = [];
            end
        end

               
        function regression_fitted(self, show, save)
            if ~isempty(self.model.fitted_estimates)
                actual = self.model.y;
                fitted = self.model.fitted_estimates;
                dates = self.complementary_information.dates;
                path = self.path;
                name = char(self.complementary_information.endogenous_variables(1));
                file_name = ['fit-' name '.png'];
                fig = gu.fit_single_variable(actual, fitted, dates, name);
                gu.show_and_save(fig, show, save, path, file_name);
            end
        end

    
        function regression_residuals(self, show, save)
            if ~isempty(self.model.residual_estimates)
                residuals = self.model.residual_estimates;
                dates = self.complementary_information.dates;
                path = self.path;
                name = char(self.complementary_information.endogenous_variables(1));
                file_name = ['residuals-' name '.png'];
                fig = gu.residual_single_variable(residuals, dates, name);
                gu.show_and_save(fig, show, save, path, file_name);
            end
        end


        function regression_forecasts(self, show, save)
            if ~isempty(self.model.forecast_estimates)
                forecasts = self.model.forecast_estimates;
                y_p = self.complementary_information.y_p;
                dates = self.complementary_information.forecast_dates;
                path = self.path;
                name = char(self.complementary_information.endogenous_variables(1));
                file_name = ['forecasts-' name '.png'];
                fig = gu.ols_forecasts_single_variable(forecasts, y_p, dates, name);
                gu.show_and_save(fig, show, save, path, file_name);
            end
        end
        
    end
    
    
end