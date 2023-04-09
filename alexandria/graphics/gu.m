classdef gu
    

    % iu stands for graphics utilities
    % A class containing static methods to handle estimation graphics 

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        function fit_single_variable(actual, fitted, dates, name, path)
            
            % fit_single_variable(actual, fitted, dates, name, path)
            % creates a plot for fit vs. actual, for a single variable
            %
            % parameters:
            % actual : array of dimension (n,1)
            %     vector of actual (observed) values
            % fitted : array of dimension (n,1)
            %     vector of fitted (in-sample predictions) values
            % dates : date array
            %     array of insample datetime entries
            % name : char
            %     name of variable being plotted as actual
            % path : char
            %     full path to folder where plot is going to be saved as image
            
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim, max_YLim] = gu.set_min_and_max([actual fitted], 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot actual and fitted
            hold on;
            plot(dates, actual, 'LineWidth', 1.3, 'color', [0.1 0.3 0.8]);
            plot(dates, fitted, 'LineWidth', 1.3, 'color', [0 0.6 0]);
            hold off;
            % set graphic limits, background and font size for ticks
            set(gca,'XLim', [dates(1) dates(end)]);
            set(gca,'YLim', [min_YLim max_YLim]);
            set(gca, 'color', [.9 .9 .9]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            ax = gca;
            new_ax = axes('Position', get(ax, 'Position'), 'FontSize', ... 
            10, 'Color', 'None', 'XColor', 'k', 'YColor', 'k', 'LineWidth', ...
            0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', ...
            'right', 'YTick', []);
            % title
            title(['Actual and fitted: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
            % save in folder as image
            full_path = fullfile(path, 'graphics', ['fit_' name '.png']);
            saveas(gcf, full_path);
            % close figure
            close(fig);
        end

        
        function residual_single_variable(residuals, dates, name, path)
            
            % residual_single_variable(residuals, dates, name, path)
            % creates a plot for the residuals, for a single variable
            %
            % parameters:
            % residuals : array of dimension (n,1)
            %     vector of residuals
            % dates : date array
            %     array of insample datetime entries
            % name : char
            %     name of variable being plotted as residuals
            % path : char
            %     full path to folder where plot is goig to be saved as image
            
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim, max_YLim] = gu.set_min_and_max(residuals, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot residuals
            hold on;
            plot(dates, residuals, 'LineWidth', 1.3, 'color', [0.1 0.3 .8]);
            plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 0.5);
            hold off;
            % set graphic limits, background and font size for ticks
            set(gca,'XLim', [dates(1) dates(end)]);
            set(gca,'YLim', [min_YLim max_YLim]);
            set(gca, 'color', [.9 .9 .9]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            ax = gca;
            new_ax = axes('Position', get(ax, 'Position'), 'FontSize', ... 
            10, 'Color', 'None', 'XColor', 'k', 'YColor', 'k', 'LineWidth', ...
            0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', ...
            'right', 'YTick', []);
            % title
            title(['Residuals: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
            % save in folder as image
            full_path = fullfile(path, 'graphics', ['residuals_' name '.png']);
            saveas(gcf, full_path);
            % close figure
            close(fig);
        end   
        
        
        function ols_forecasts_single_variable(forecasts, y_p, dates, name, path)
            
            % ols_forecasts_single_variable(forecasts, dates, name, path)
            % creates a plot for the forecasts of an OLS model
            %
            % parameters:
            % forecasts : array
            %     array of forecast estimates
            % y_p : array
            %     array (possibly empty) of actual values for the forecasts
            % dates : date array
            %     array of forecast datetime entries
            % name : char
            %     name of variable being plotted as forecast
            % path : char
            %     full path to folder where plot is goig to be saved as image
            
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim, max_YLim] = gu.set_min_and_max([forecasts y_p], 0.15, 0.15);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot forecasts, with patched credibility intervals
            hold on;
            x_patch = [dates;flipud(dates)];
            y_patch = [forecasts(:,1);flipud(forecasts(:,3))];
            forecast_patch = fill(x_patch, y_patch, 'r');
            set(forecast_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', ...
                'none', 'facealpha', 0.4);
            plot(dates, forecasts(:,2), 'LineWidth', 1.6, 'color', [0 0.6 0]);
            if y_p
                plot(dates, y_p, 'LineWidth', 1.6, 'color', [0.1 0.3 0.8]);
            end
            hold off;
            % set graphic limits, background and font size for ticks
            set(gca,'XLim', [dates(1) dates(end)]);
            set(gca,'YLim', [min_YLim max_YLim]);
            set(gca, 'color', [.92 .92 .92]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            ax = gca;
            new_ax = axes('Position', get(ax, 'Position'), 'FontSize', ... 
            10, 'Color', 'None', 'XColor', 'k', 'YColor', 'k', 'LineWidth', ...
            0.5, 'XAxisLocation', 'top', 'XTick', [], 'YAxisLocation', ...
            'right', 'YTick', []);
            % title
            title(['Forecasts: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
            % save in folder as image
            full_path = fullfile(path, 'graphics', ['forecasts_' name '.png']);
            saveas(gcf, full_path);
            % close figure
            close(fig);
        end           
        
        
        function [min_YLim, max_YLim] = set_min_and_max(data, min_space, max_space)
            
            % [min_YLim, max_YLim] = set_min_and_max(data)
            % returns the min and max Y values for a graphic, given data
            %
            % parameters:
            % data : array
            %     matrix of plotted values
            % min_space : float
            %     scale to define lower space
            % max_space : float
            %     scale to define upper space
            % 
            % returns:
            % min_YLim : scalar
            %     min Y value for the plot
            % max_YLim : scalar
            %     max Y value for the plot
            
            % get min, max, and compute window width
            min_value = min(data, [], 'all');
            max_value = max(data, [], 'all');
            width = max_value - min_value;
            min_YLim = min_value - min_space * width;
            max_YLim = max_value + max_space * width;
        end
  

    end
        
        
end
    
    
    