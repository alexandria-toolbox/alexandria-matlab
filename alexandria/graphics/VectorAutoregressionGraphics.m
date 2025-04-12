classdef VectorAutoregressionGraphics < handle
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)   


        function [self] = VectorAutoregressionGraphics()    
        end


    end
    
    
    methods (Access = protected, Hidden = true)


        function complete_var_information(self)  
            % endogenous and exogenous variables
            if ~isfield(self.complementary_information, 'endogenous_variables')
                n_endo = size(self.model.endogenous,2);
                self.complementary_information.endogenous_variables = "y"+(1:n_endo);
            end       
            if ~isfield(self.complementary_information, 'exogenous_variables')
                if ~isempty(self.model.exogenous)
                    n_exo = size(self.model.exogenous,2);
                    self.complementary_information.exogenous_variables = "x"+(1:n_exo);
                else
                    self.complementary_information.exogenous_variables = [];
                end
            end
            % sample dates
            if ~isfield(self.complementary_information, 'dates')
                T = self.model.T;
                p = self.model.p;
                self.complementary_information.dates = (-p+1:T)';
            end
            % forecast_dates
            if ~isempty(self.model.forecast_estimates) && ~isfield(self.complementary_information, 'forecast_dates')
                T = self.model.T;
                forecast_periods = size(self.model.forecast_estimates,1);
                forecast_dates = (T+1:T+forecast_periods)';
                self.complementary_information.forecast_dates = forecast_dates;
            end
            % conditional forecast_dates
            if ~isempty(self.model.conditional_forecast_estimates) && ~isfield(self.complementary_information, 'conditional_forecast_dates')
                T = self.model.T;
                forecast_periods = size(self.model.conditional_forecast_estimates,1);
                forecast_dates = (T+1:T+forecast_periods)';
                self.complementary_information.conditional_forecast_dates = forecast_dates;
            end  
            % actual
            if ~isfield(self.complementary_information, 'Y_p')
                self.complementary_information.Y_p = [];
            end
            % structural shocks
            if ~isfield(self.complementary_information, 'shocks')
                n = self.model.n;
                self.complementary_information.shocks = "shock"+(1:n);
            end
        end


        function var_fitted(self, show, save)
            if ~isempty(self.model.fitted_estimates)
                % recover graphics elements
                n = self.model.n;
                p = self.model.p;
                actual = self.model.Y;
                fitted = permute(self.model.fitted_estimates, [1 3 2]);
                dates = self.complementary_information.dates(p+1:end);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;
                % produce individual graphs
                for i=1:n
                    file_name = ['fit-' char(endogenous(i)) '.png'];
                    fig = gu.var_fit_single_variable(actual(:,i), fitted(:,:,i), dates, endogenous(i));
                    gu.show_and_save(fig, show, save, path, file_name);
                end
                % joint graph
                fig = gu.var_fit_all(actual, fitted, dates, endogenous, n);
                gu.show_and_save(fig, show, save, path, 'fit-all.png') ;
            end
        end


        function var_residuals(self, show, save)
            if ~isempty(self.model.residual_estimates)
                % recover graphics elements
                n = self.model.n;
                p = self.model.p;
                residuals = permute(self.model.residual_estimates, [1 3 2]);
                dates = self.complementary_information.dates(p+1:end);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;          
                % produce individual graphs
                for i=1:n
                    file_name = ['residuals-' char(endogenous(i)) '.png'];
                    fig = gu.var_residual_single_variable(residuals(:,:,i), dates, endogenous(i));
                    gu.show_and_save(fig, show, save, path, file_name);
                end
                % joint graph
                fig = gu.var_residual_all(residuals, dates, endogenous, n);
                gu.show_and_save(fig, show, save, path, 'residuals-all.png');
            end
        end


        function var_shocks(self, show, save)  
            if ~isempty(self.model.structural_shock_estimates)
                % recover graphics elements
                n = self.model.n;
                p = self.model.p;
                shocks = permute(self.model.structural_shock_estimates, [1 3 2]);
                dates = self.complementary_information.dates(p+1:end);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;            
                % produce individual graphs
                for i=1:n
                    file_name = ['shocks-' char(endogenous(i)) '.png'];
                    fig = gu.var_shocks_single_variable(shocks(:,:,i), dates, endogenous(i));
                    gu.show_and_save(fig, show, save, path, file_name);
                end
                % joint graph
                fig = gu.var_shocks_all(shocks, dates, endogenous, n);
                gu.show_and_save(fig, show, save, path, 'shocks-all.png');
            end
        end


        function var_steady_state(self, show, save)  
            if ~isempty(self.model.steady_state_estimates)
                % recover graphics elements
                n = self.model.n;
                p = self.model.p;
                actual = self.model.Y;
                steady_state = permute(self.model.steady_state_estimates, [1 3 2]);
                dates = self.complementary_information.dates(p+1:end);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;             
                % produce individual graphs
                for i=1:n
                    file_name = ['steady_state-' char(endogenous(i)) '.png'];
                    fig = gu.var_steady_state_single_variable(actual(:,i),steady_state(:,:,i), dates, endogenous(i));
                    gu.show_and_save(fig, show, save, path, file_name);
                end
                % joint graph
                fig = gu.var_steady_state_all(actual, steady_state, dates, endogenous, n);
                gu.show_and_save(fig, show, save, path, 'steady_state-all.png');
            end
        end


        function var_forecasts(self, show, save)
            if ~isempty(self.model.forecast_estimates)
                % recover graphics elements
                n = self.model.n;
                p = self.model.p;
                actual = self.model.Y;
                forecasts = permute(self.model.forecast_estimates, [1 3 2]);
                Y_p = self.complementary_information.Y_p;
                Y_p_i = [];
                dates = self.complementary_information.dates(p+1:end);
                forecast_dates = self.complementary_information.forecast_dates;
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;
                % produce individual graphs
                for i=1:n
                    file_name = ['forecasts-' char(endogenous(i)) '.png'];
                    if ~isempty(Y_p)
                        Y_p_i = Y_p(:,i);
                    end
                    fig = gu.var_forecasts_single_variable(actual(:,i),forecasts(:,:,i), Y_p_i, dates, forecast_dates, endogenous(i));
                    gu.show_and_save(fig, show, save, path, file_name);   
                end
                % joint graph
                fig = gu.var_forecasts_all(actual, forecasts, Y_p, dates, forecast_dates, endogenous, n);
                gu.show_and_save(fig, show, save, path, 'forecasts-all.png');
            end
        end


        function var_conditional_forecasts(self, show, save)
            if isprop(self.model,'conditional_forecast_estimates') && ...
               ~isempty(self.model.conditional_forecast_estimates)
                % recover graphics elements
                n = self.model.n;
                p = self.model.p;
                actual = self.model.Y;
                forecasts = permute(self.model.conditional_forecast_estimates, [1 3 2]);
                Y_p = self.complementary_information.Y_p;
                Y_p_i = [];
                dates = self.complementary_information.dates(p+1:end);
                forecast_dates = self.complementary_information.conditional_forecast_dates;
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;
                % produce individual graphs
                for i=1:n
                    file_name = ['conditional_forecasts-' char(endogenous(i)) '.png'];
                    if ~isempty(Y_p)
                        Y_p_i = Y_p(:,i);
                    end
                    fig = gu.var_conditional_forecasts_single_variable(actual(:,i),...
                          forecasts(:,:,i), Y_p_i, dates, forecast_dates, endogenous(i));
                    gu.show_and_save(fig, show, save, path, file_name);   
                end
                % joint graph
                fig = gu.var_conditional_forecasts_all(actual, forecasts, Y_p, ...
                      dates, forecast_dates, endogenous, n);
                gu.show_and_save(fig, show, save, path, 'conditional_forecasts-all.png');
            end
        end


        function var_irf(self, show, save)
            if ~isempty(self.model.irf_estimates)  
                % recover graphics elements
                n = self.model.n;
                irf = permute(self.model.irf_estimates, [3 4 1 2]);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;
                shocks = self.complementary_information.shocks;
                % produce individual graphs
                for i=1:n
                    for j=1:n
                        file_name = ['irf-' char(endogenous(i)) '@'  char(shocks(j)) '.png'];
                        fig = gu.var_irf_single_variable(irf(:,:,i,j), char(endogenous(i)), char(shocks(j)));
                        gu.show_and_save(fig, show, save, path, file_name);
                    end
                end
            end
            if ~isempty(self.model.exo_irf_estimates)
                % recover graphics elements
                n = self.model.n;
                n_exo = size(self.model.exogenous, 2);
                exo_irf = permute(self.model.exo_irf_estimates, [3 4 1 2]);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;
                shocks = self.complementary_information.exogenous_variables;
                % produce individual graphs
                for i=1:n
                    for j=1:n_exo
                        file_name = ['irf-' char(endogenous(i)) '@' char(shocks(j)) '.png'];
                        fig = gu.var_irf_single_variable(exo_irf(:,:,i,j), char(endogenous(i)), char(shocks(j)));
                        gu.show_and_save(fig, show, save, path, file_name);
                    end
                end
            end
            if ~isempty(self.model.irf_estimates)
                % recover graphics elements
                n_endo = self.model.n;
                n_shocks = n_endo;
                irf = permute(self.model.irf_estimates, [3 4 1 2]);
                variables = self.complementary_information.endogenous_variables;
                shocks = self.complementary_information.shocks;
                if ~isempty(self.model.exo_irf_estimates)
                    n_exo = size(self.model.exogenous, 2);
                    exo_irf = permute(self.model.exo_irf_estimates, [3 4 1 2]);
                    exogenous = self.complementary_information.exogenous_variables;
                    n_shocks = n_shocks + n_exo;
                    irf = cat(4, irf, exo_irf);
                    shocks = [shocks exogenous];
                end
                % joint graph
                fig = gu.var_irf_all(irf, variables, shocks, n_endo, n_shocks);
                gu.show_and_save(fig, show, save, path, 'irf-all.png');
            end
        end


        function var_fevd(self, show, save)
            if ~isempty(self.model.fevd_estimates)
                % recover graphics elements
                n = self.model.n;
                fevd = permute(self.model.fevd_estimates, [3 4 1 2]);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;
                shocks = self.complementary_information.shocks;
                % produce individual graphs
                for i=1:n
                    for j=1:n
                        file_name = ['fevd-' char(endogenous(i)) '@' char(shocks(j)) '.png'];
                        fig = gu.var_fevd_single_variable(fevd(:,:,i,j), char(endogenous(i)), char(shocks(j)));
                        gu.show_and_save(fig, show, save, path, file_name);
                    end
                end
                % partial joint graph
                fevd = permute(self.model.fevd_estimates, [3 2 1 4]);
                for i=1:n
                    file_name = ['fevd-' char(endogenous(i)) '@all.png'];
                    fig = gu.var_fevd_joint(fevd(:,:,i,1), char(endogenous(i)), shocks, n);
                    gu.show_and_save(fig, show, save, path, file_name);
                end
                % joint graph
                fevd = permute(self.model.fevd_estimates, [3 2 1 4]);
                fig = gu.var_fevd_all(fevd, endogenous, shocks, n);
                gu.show_and_save(fig, show, save, path, 'fevd-all.png');
            end
        end


        function var_hd(self, show, save)
            if ~isempty(self.model.hd_estimates)
                % recover graphics elements
                n = self.model.n;
                p = self.model.p;
                T = self.model.T;
                dates = self.complementary_information.dates(p+1:end);
                hd = permute(self.model.hd_estimates, [3 4 1 2]);
                path = self.path;
                endogenous = self.complementary_information.endogenous_variables;
                shocks = self.complementary_information.shocks;
                % produce individual graphs
                for i=1:n
                    for j=1:n
                        file_name = ['hd-' char(endogenous(i)) '@' char(shocks(j)) '.png'];
                        fig = gu.var_hd_single_variable(hd(:,:,i,j), char(endogenous(i)), char(shocks(j)), dates);
                        gu.show_and_save(fig, show, save, path, file_name);
                    end
                end
                % partial joint graph
                hd = permute(self.model.hd_estimates, [3 2 1 4]);
                for i=1:n
                    file_name = ['hd-' char(endogenous(i)) '@all.png'];
                    fig = gu.var_hd_joint(hd(:,:,i,1), char(endogenous(i)), shocks, dates, n, T);
                    gu.show_and_save(fig, show, save, path, file_name);
                end
                % joint graph
                hd = permute(self.model.hd_estimates, [3 2 1 4]);
                fig = gu.var_hd_all(hd, endogenous, shocks, dates, n, T);
                gu.show_and_save(fig, show, save, path, 'hd-all.png');
            end
        end

    end
    
end