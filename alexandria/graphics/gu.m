classdef gu
    

    % gu stands for graphics utilities
    % a class containing static methods to handle estimation graphics 

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        

        function [fig] = fit_single_variable(actual, fitted, dates, name)

            % fit_single_variable(actual, fitted, dates, name)
            % produces fitted figure for regression model, single variable
            %
            % parameters:
            % actual : array of dimension (n,1)
            %     actual sample values
            % fitted : array of dimension (n,1)
            %     fitted values
            % dates : date array
            %     array of in-sample dates
            % name : char
            %     name of variable for which figure is produced
            % 
            % returns:
            % fig : matlab figure
            %     fitted figure

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max([actual fitted], 0.07, 0.1);
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
        end

        
        function [fig] = residual_single_variable(residuals, dates, name)
            
            % residual_single_variable(residuals, dates, name)
            % produces residual figure for regression model, single variable
            % 
            % residuals: matrix of size (n,1)
            %     residual values
            % dates: datetime array of size (n,1)
            %     index of in-sample dates
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     residual figure
            
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max(residuals, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot residuals
            hold on;
            plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
            plot(dates, residuals, 'LineWidth', 1.3, 'color', [0.1 0.3 .8]);
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
        end   
        

        function [fig] = ols_forecasts_single_variable(forecasts, y_p, dates, name)
            
            % ols_forecasts_single_variable(forecasts, y_p, dates, name)
            % produces forecast figure for regression model, single variable
            % 
            % forecasts: matrix of size (n_forecast,3)
            %     forecast values, median, lower and upper bound
            % y_p: ndarray of size (n_forecast,) or empty array
            %     actual values for forecast evaluation
            % dates: datetime array of size (n,1)
            %     index of in-sample dates
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     forecast figure

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max([forecasts y_p], 0.15, 0.15);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot forecasts, with patched credibility intervals
            hold on;
            x_patch = [dates;flipud(dates)];
            y_patch = [forecasts(:,2);flipud(forecasts(:,3))];
            forecast_patch = fill(x_patch, y_patch, 'r');
            set(forecast_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', ...
                'none', 'facealpha', 0.4);
            if min_YLim < 0 && max_YLim > 0
                plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
            end
            plot(dates, forecasts(:,1), 'LineWidth', 1.6, 'color', [0 0.6 0]);
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
        end 


        function [fig] = var_fit_single_variable(actual, fitted, dates, name)

            % var_fit_single_variable(actual, fitted, dates, name)
            % produces fitted figure for var model, single variable
            % 
            % parameters:     
            % actual: matrix of size (T,1)
            %     actual sample values
            % fitted: matrix of size (T,3)
            %     fitted values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max([actual fitted], 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot actual and fitted
            hold on;
            x_patch = [dates;flipud(dates)];
            y_patch = [fitted(:,2);flipud(fitted(:,3))];
            fitted_patch = fill(x_patch, y_patch, 'g');
            set(fitted_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot(dates, fitted(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            plot(dates, actual, 'LineWidth', 1.3, 'color', [0.1 0.3 0.8]);
            plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
            % title
            title(['Actual and fitted: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_fit_all(actual, fitted, dates, endogenous, n)

            % var_fit_all(actual, fitted, dates, endogenous, n)
            % produces fitted figure for var model, all variables
            % 
            % parameters:     
            % actual: matrix of size (T,n)
            %     actual sample values
            % fitted: matrix of size (T,n,3)
            %     fitted values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % endogenous: str array
            %     list of endogenous variables for which figure is produced
            % n: int
            %     number of endogenous variables            
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            x_patch = [dates;flipud(dates)];
            for i = 1:n
                axes('Units', 'normalized', 'Position', positions(i,:));
                % get min and max for subplot
                [min_YLim max_YLim] = gu.set_min_and_max([actual(:,i) fitted(:,:,i)], 0.07, 0.1);
                % plot actual, and fitted with patched credibility intervals
                hold on;
                y_patch = [fitted(:,2,i);flipud(fitted(:,3,i))];
                fitted_patch = fill(x_patch, y_patch, 'g');
                set(fitted_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                plot(dates, fitted(:,1,i), 'LineWidth', 1.3, 'color', [0 0.5 0]);
                plot(dates, actual(:,i), 'LineWidth', 1.3, 'color', [0.1 0.3 0.8]);
                plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
                % title
                name = char(endogenous(i));
                title(['Actual and fitted: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [fig] = var_residual_single_variable(residuals, dates, name)

            % var_fit_single_variable(actual, fitted, dates, name)
            % produces fitted figure for var model, single variable
            % 
            % parameters:     
            % residuals: matrix of size (T,3)
            %     residual values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max(residuals, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot residuals with patched credibility intervals
            hold on;
            x_patch = [dates;flipud(dates)];
            y_patch = [residuals(:,2);flipud(residuals(:,3))];
            fitted_patch = fill(x_patch, y_patch, 'g');
            set(fitted_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
            plot(dates, residuals(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
            % title
            title(['Residuals: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_residual_all(residuals, dates, endogenous, n)

            % var_residual_all(residuals, dates, endogenous, n)
            % produces residual figure for var model, all variables
            % 
            % parameters:     
            % residuals: matrix of size (T,n,3)
            %     residual values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % endogenous: str array
            %     list of endogenous variables for which figure is produced
            % n: int
            %     number of endogenous variables   
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            x_patch = [dates;flipud(dates)];
            for i = 1:n
                axes('Units', 'normalized', 'Position', positions(i,:));
                % get min and max for subplot
                [min_YLim max_YLim] = gu.set_min_and_max([residuals(:,:,i)], 0.07, 0.1);
                % plot residuals with patched credibility intervals
                hold on;
                y_patch = [residuals(:,2,i);flipud(residuals(:,3,i))];
                residual_patch = fill(x_patch, y_patch, 'g');
                set(residual_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
                plot(dates, residuals(:,1,i), 'LineWidth', 1.3, 'color', [0 0.5 0]);
                plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
                % title
                name = char(endogenous(i));
                title(['Residuals: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [fig] = var_shocks_single_variable(shocks, dates, name)

            % var_shocks_single_variable(shocks, dates, name)
            % produces structural shocks figure for var model, single variable
            % 
            % parameters:     
            % shocks: matrix of size (T,3)
            %     shock values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max(shocks, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot shocks with patched credibility intervals
            hold on;
            x_patch = [dates;flipud(dates)];
            y_patch = [shocks(:,2);flipud(shocks(:,3))];
            shock_patch = fill(x_patch, y_patch, 'g');
            set(shock_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
            plot(dates, shocks(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
            % title
            title(['Shocks: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_shocks_all(shocks, dates, endogenous, n)

            % var_shocks_all(shocks, dates, endogenous, n)
            % produces shock figure for var model, all variables
            % 
            % parameters:     
            % shocks: matrix of size (T,n,3)
            %     shock values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % endogenous: str array
            %     list of endogenous variables for which figure is produced
            % n: int
            %     number of endogenous variables   
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            x_patch = [dates;flipud(dates)];
            for i = 1:n
                axes('Units', 'normalized', 'Position', positions(i,:));
                % get min and max for subplot
                [min_YLim max_YLim] = gu.set_min_and_max([shocks(:,:,i)], 0.07, 0.1);
                % plot shocks with patched credibility intervals
                hold on;
                y_patch = [shocks(:,2,i);flipud(shocks(:,3,i))];
                shock_patch = fill(x_patch, y_patch, 'g');
                set(shock_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
                plot(dates, shocks(:,1,i), 'LineWidth', 1.3, 'color', [0 0.5 0]);
                plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
                % title
                name = char(endogenous(i));
                title(['Shocks: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [fig] = var_steady_state_single_variable(actual, steady_state, dates, name)

            % var_steady_state_single_variable(actual, steady_state, dates, name)
            % produces steady-state figure for var model, single variable
            % 
            % parameters: 
            % actual: matrix of size (T,1)
            %     actual sample values
            % steady_state: matrix of size (T,3)
            %     steady-state values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max([actual steady_state], 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot shocks with patched credibility intervals
            hold on;
            x_patch = [dates;flipud(dates)];
            y_patch = [steady_state(:,2);flipud(steady_state(:,3))];
            steady_patch = fill(x_patch, y_patch, 'g');
            set(steady_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot(dates, steady_state(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            plot(dates, actual, 'LineWidth', 1.3, 'color', [0.1 0.3 0.8]);
            plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
            % title
            title(['Steady-state: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_steady_state_all(actual, steady_state, dates, endogenous, n)

            % var_steady_state_all(actual, steady_state, dates, endogenous, n)
            % produces steady-state figure for var model, all variables
            % 
            % parameters:   
            % actual: matrix of size (T,n)
            %     actual sample values
            % steady_state: matrix of size (T,n,3)
            %     steady_state values, median, lower and upper bounds
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % endogenous: str array
            %     list of endogenous variables for which figure is produced
            % n: int
            %     number of endogenous variables   
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            x_patch = [dates;flipud(dates)];
            for i = 1:n
                axes('Units', 'normalized', 'Position', positions(i,:));
                % get min and max for subplot
                [min_YLim max_YLim] = gu.set_min_and_max([actual(:,i) steady_state(:,:,i)], 0.07, 0.1);
                % plot actual, and steady-state with patched credibility intervals
                hold on;
                y_patch = [steady_state(:,2,i);flipud(steady_state(:,3,i))];
                steady_patch = fill(x_patch, y_patch, 'g');
                set(steady_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                plot(dates, steady_state(:,1,i), 'LineWidth', 1.3, 'color', [0 0.5 0]);
                plot(dates, actual(:,i), 'LineWidth', 1.3, 'color', [0.1 0.3 0.8]);
                plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
                % title
                name = char(endogenous(i));
                title(['Steady-state: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [fig] = var_forecasts_single_variable(actual, forecasts, Y_p, dates, forecast_dates, name)

            % var_forecasts_single_variable(actual, forecasts, Y_p, dates, forecast_dates, name)
            % produces forecast figure for var model, single variable
            % 
            % parameters: 
            % actual: matrix of size (T,1)
            %     actual sample values
            % forecasts: matrix of size (f_periods,3)
            %     forecast values, median, lower and upper bounds
            % Y_p: matrix of size (f_periods,1)
            %     actual out-of-sample values  
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % forecast_dates: datetime array of size (f_periods)
            %     index of forecast dates            
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % periods to plot
            T = max(10, 2 * size(forecasts,1));
            sample_data = actual(end-T+1:end);
            sample_dates = dates(end-T+1:end);
            if isempty(Y_p)
                plot_data = [repmat(sample_data,[1 3]);forecasts];
            else
                plot_data = [repmat(sample_data, [1 4]);[forecasts Y_p]];
            end
            plot_dates = [sample_dates;forecast_dates];
            prediction_data = plot_data(T:end,:);
            prediction_dates = plot_dates(T:end);
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max(plot_data, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot actual, and forecasts with patched credibility intervals
            hold on;
            x_patch = [prediction_dates;flipud(prediction_dates)];
            y_patch = [prediction_data(:,2);flipud(prediction_data(:,3))];
            forecast_patch = fill(x_patch, y_patch, 'g');
            set(forecast_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot(prediction_dates, prediction_data(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            if ~isempty(Y_p)
                plot(prediction_dates, prediction_data(:,4), 'LineWidth', 1.3, 'LineStyle', '--', 'color', [0.1 0.3 0.8]);
            end
            plot(sample_dates, sample_data, 'Linewidth', 1.5, 'color', [0.1, 0.3, 0.8]);
            plot([plot_dates(1) plot_dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([plot_dates(end) plot_dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            hold off;
            % set graphic limits, background and font size for ticks
            set(gca,'XLim', [plot_dates(1) plot_dates(end)]);
            set(gca,'YLim', [min_YLim max_YLim]);
            set(gca, 'color', [.9 .9 .9]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            % title
            title(['Forecasts: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_forecasts_all(actual, forecasts, Y_p, dates, forecast_dates, endogenous, n)

            % var_forecasts_all(actual, forecasts, Y_p, dates, forecast_dates, endogenous, n)
            % produces forecast figure for var model, all variables
            % 
            % parameters:     
            % shocks: matrix of size (T,n,3)
            %     shock values, median, lower and upper bounds
            % actual: matrix of size (T,n)
            %     actual sample values
            % forecasts: matrix of size (f_periods,n,3)
            %     forecast values, median, lower and upper bounds
            % Y_p: matrix of size (f_periods,n)
            %     actual out-of-sample values
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % forecast_dates: datetime array of size (f_periods)
            %     index of forecast dates                     
            % endogenous: str array
            %     list of endogenous variables for which figure is produced
            % n: int
            %     number of endogenous variables   
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % periods to plot
            T = max(10, 2 * size(forecasts,1));
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            for i = 1:n
                sample_data = actual(end-T+1:end,i);
                sample_dates = dates(end-T+1:end);
                if isempty(Y_p)
                    plot_data = [repmat(sample_data,[1 3]);forecasts(:,:,i)];
                else
                    plot_data = [repmat(sample_data, [1 4]);[forecasts(:,:,i) Y_p(:,i)]];
                end                
                plot_dates = [sample_dates;forecast_dates];
                prediction_data = plot_data(T:end,:);
                prediction_dates = plot_dates(T:end);
                axes('Units', 'normalized', 'Position', positions(i,:));
                % get min and max for subplot
                [min_YLim max_YLim] = gu.set_min_and_max(plot_data, 0.07, 0.1);
                % plot actual, and steady-state with patched credibility intervals
                hold on;
                x_patch = [prediction_dates;flipud(prediction_dates)];
                y_patch = [prediction_data(:,2);flipud(prediction_data(:,3))];
                forecast_patch = fill(x_patch, y_patch, 'g');
                set(forecast_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                plot(prediction_dates, prediction_data(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
                if ~isempty(Y_p)
                    plot(prediction_dates, prediction_data(:,4), 'LineWidth', 1.3, 'LineStyle', '--', 'color', [0.1 0.3 0.8]);
                end 
                plot(sample_dates, sample_data, 'Linewidth', 1.5, 'color', [0.1, 0.3, 0.8]);
                plot([plot_dates(1) plot_dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([plot_dates(end) plot_dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                hold off;
                % set graphic limits, background and font size for ticks
                set(gca,'XLim', [plot_dates(1) plot_dates(end)]);
                set(gca,'YLim', [min_YLim max_YLim]);
                set(gca, 'color', [.9 .9 .9]);
                set(gca, 'FontSize',13);
                % set figure and plot background color
                set(0,'defaultfigurecolor',[1 1 1]);
                grid on;
                set(gca, 'GridColor', [.3 .3 .3]);
                % create top and right axes
                box off;
                % title
                name = char(endogenous(i));
                title(['Forecasts: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [fig] = var_conditional_forecasts_single_variable(actual, forecasts, Y_p, dates, forecast_dates, name)

            % var_conditional_forecasts_single_variable(actual, forecasts, Y_p, dates, forecast_dates, name)
            % produces conditional forecast figure for var model, single variable
            % 
            % parameters: 
            % actual: matrix of size (T,1)
            %     actual sample values
            % forecasts: matrix of size (f_periods,3)
            %     forecast values, median, lower and upper bounds
            % Y_p: matrix of size (f_periods,1)
            %     actual out-of-sample values  
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % forecast_dates: datetime array of size (f_periods)
            %     index of forecast dates            
            % name: char
            %     name of variable for which figure is produced
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % periods to plot
            T = max(10, 2 * size(forecasts,1));
            sample_data = actual(end-T+1:end);
            sample_dates = dates(end-T+1:end);
            if isempty(Y_p)
                plot_data = [repmat(sample_data,[1 3]);forecasts];
            else
                plot_data = [repmat(sample_data, [1 4]);[forecasts Y_p]];
            end
            plot_dates = [sample_dates;forecast_dates];
            prediction_data = plot_data(T:end,:);
            prediction_dates = plot_dates(T:end);
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max(plot_data, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot actual, and forecasts with patched credibility intervals
            hold on;
            x_patch = [prediction_dates;flipud(prediction_dates)];
            y_patch = [prediction_data(:,2);flipud(prediction_data(:,3))];
            forecast_patch = fill(x_patch, y_patch, 'g');
            set(forecast_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot(prediction_dates, prediction_data(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            if ~isempty(Y_p)
                plot(prediction_dates, prediction_data(:,4), 'LineWidth', 1.3, 'LineStyle', '--', 'color', [0.1 0.3 0.8]);
            end
            plot(sample_dates, sample_data, 'Linewidth', 1.5, 'color', [0.1, 0.3, 0.8]);
            plot([plot_dates(1) plot_dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([plot_dates(end) plot_dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            hold off;
            % set graphic limits, background and font size for ticks
            set(gca,'XLim', [plot_dates(1) plot_dates(end)]);
            set(gca,'YLim', [min_YLim max_YLim]);
            set(gca, 'color', [.9 .9 .9]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            % title
            title(['Conditional forecasts: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_conditional_forecasts_all(actual, forecasts, Y_p, dates, forecast_dates, endogenous, n)

            % var_conditional_forecasts_all(actual, forecasts, Y_p, dates, forecast_dates, endogenous, n)
            % produces conditional forecast figure for var model, all variables
            % 
            % parameters:     
            % shocks: matrix of size (T,n,3)
            %     shock values, median, lower and upper bounds
            % actual: matrix of size (T,n)
            %     actual sample values
            % forecasts: matrix of size (f_periods,n,3)
            %     forecast values, median, lower and upper bounds
            % Y_p: matrix of size (f_periods,n)
            %     actual out-of-sample values
            % dates: datetime array of size (T)
            %     index of in-sample dates
            % forecast_dates: datetime array of size (f_periods)
            %     index of forecast dates                     
            % endogenous: str array
            %     list of endogenous variables for which figure is produced
            % n: int
            %     number of endogenous variables   
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % periods to plot
            T = max(10, 2 * size(forecasts,1));
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            for i = 1:n
                sample_data = actual(end-T+1:end,i);
                sample_dates = dates(end-T+1:end);
                if isempty(Y_p)
                    plot_data = [repmat(sample_data,[1 3]);forecasts(:,:,i)];
                else
                    plot_data = [repmat(sample_data, [1 4]);[forecasts(:,:,i) Y_p(:,i)]];
                end                
                plot_dates = [sample_dates;forecast_dates];
                prediction_data = plot_data(T:end,:);
                prediction_dates = plot_dates(T:end);
                axes('Units', 'normalized', 'Position', positions(i,:));
                % get min and max for subplot
                [min_YLim max_YLim] = gu.set_min_and_max(plot_data, 0.07, 0.1);
                % plot actual, and steady-state with patched credibility intervals
                hold on;
                x_patch = [prediction_dates;flipud(prediction_dates)];
                y_patch = [prediction_data(:,2);flipud(prediction_data(:,3))];
                forecast_patch = fill(x_patch, y_patch, 'g');
                set(forecast_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                plot(prediction_dates, prediction_data(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
                if ~isempty(Y_p)
                    plot(prediction_dates, prediction_data(:,4), 'LineWidth', 1.3, 'LineStyle', '--', 'color', [0.1 0.3 0.8]);
                end 
                plot(sample_dates, sample_data, 'Linewidth', 1.5, 'color', [0.1, 0.3, 0.8]);
                plot([plot_dates(1) plot_dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([plot_dates(end) plot_dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                hold off;
                % set graphic limits, background and font size for ticks
                set(gca,'XLim', [plot_dates(1) plot_dates(end)]);
                set(gca,'YLim', [min_YLim max_YLim]);
                set(gca, 'color', [.9 .9 .9]);
                set(gca, 'FontSize',13);
                % set figure and plot background color
                set(0,'defaultfigurecolor',[1 1 1]);
                grid on;
                set(gca, 'GridColor', [.3 .3 .3]);
                % create top and right axes
                box off;
                % title
                name = char(endogenous(i));
                title(['Conditional forecasts: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [fig] = var_irf_single_variable(irf, name, shock)

            % var_irf_single_variable(irf, variable, shock)
            % produces IRF figure for var model, single variable
            % 
            % parameters:     
            % irf: matrix of size (irf_periods,n,3)
            %     IRF values, median, lower and upper bounds
            % name: char
            %     name of variable for which figure is produced
            % shock: char
            %     name of shock for which figure is produced        
            %     
            % returns:
            % fig: matlab figure
            %     IRF figure

            % periods to plot
            irf_periods = (1:size(irf,1))';
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max(irf, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot irf with patched credibility intervals
            hold on;
            x_patch = [irf_periods;flipud(irf_periods)];
            y_patch = [irf(:,2);flipud(irf(:,3))];
            fitted_patch = fill(x_patch, y_patch, 'g');
            set(fitted_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot([irf_periods(1) irf_periods(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
            plot(irf_periods, irf(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            plot([irf_periods(1) irf_periods(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([irf_periods(end) irf_periods(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            hold off;
            % set graphic limits, background and font size for ticks
            set(gca,'XLim', [irf_periods(1) irf_periods(end)]);
            set(gca,'YLim', [min_YLim max_YLim]);
            set(gca, 'color', [.9 .9 .9]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            % title
            title(['IRF: ' name '\_' shock], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_irf_all(irf, variables, shocks, n_endo, n_shocks)
            
            % var_irf_all(irf, variables, shocks, n_endo, n_shocks)
            % produces IRF figure for var model, all variables
            % 
            % parameters:     
            % irf: matrix of size (irf_periods,3,n,n)
            %     IRF values, median, lower and upper bounds
            % variables: str array
            %     list of endogenous variables for which figure is produced
            % shocks: str array
            %     list of structural shocks for which figure is produced    
            % n_endo: int
            %     number of endogenous variables
            % n_exo: int
            %     number of structural shocks and exogenous variables        
            %     
            % returns:
            % fig: matlab figure
            %     IRF figure

            % periods to plot
            irf_periods = (1:size(irf,1))';
            % create figure
            fig = figure('Position', [100 100 660*n_shocks 470*n_shocks], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n_endo * n_shocks, n_shocks, n_shocks, 0.07, 0.05, 0.04, 0.05);
            x_patch = [irf_periods;flipud(irf_periods)];
            % loop over variables and shocks
            k = 0;            
            for i = 1:n_endo
                for j = 1:n_shocks
                    k = k + 1;
                    axes('Units', 'normalized', 'Position', positions(k,:));
                    % get min and max for subplot
                    [min_YLim max_YLim] = gu.set_min_and_max([irf(:,:,i,j)], 0.07, 0.1);
                    % plot IRF with patched credibility intervals
                    hold on;
                    y_patch = [irf(:,2,i,j);flipud(irf(:,3,i,j))];
                    irf_patch = fill(x_patch, y_patch, 'g');
                    set(irf_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                    plot([irf_periods(1) irf_periods(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
                    plot(irf_periods, irf(:,1,i,j), 'LineWidth', 1.3, 'color', [0 0.5 0]);
                    plot([irf_periods(1) irf_periods(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                    plot([irf_periods(end) irf_periods(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                    hold off;
                    % set graphic limits, background and font size for ticks
                    set(gca,'XLim', [irf_periods(1) irf_periods(end)]);
                    set(gca,'YLim', [min_YLim max_YLim]);
                    set(gca, 'color', [.9 .9 .9]);
                    set(gca, 'FontSize',13);
                    % set figure and plot background color
                    set(0,'defaultfigurecolor',[1 1 1]);
                    grid on;
                    set(gca, 'GridColor', [.3 .3 .3]);
                    % create top and right axes
                    box off;
                    % title
                    name = char(variables(i));
                    shock = char(shocks(j));
                    title(['IRF: ' name '\_' shock], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                    % command to preserve grey background color when saving as image
                    fig.InvertHardcopy = 'off';
                end
            end
        end


        function [fig] = var_fevd_single_variable(fevd, name, shock)

            % var_fevd_single_variable(fevd, name, shock)
            % produces FEVD figure for var model, single variable
            % 
            % parameters:     
            % fevd: matrix of size (fevd_periods,3)
            %     fevd values, median, lower and upper bounds
            % name: char
            %     name of variable for which figure is produced
            % shock: char
            %     name of shock for which figure is produced        
            %     
            % returns:
            % fig: matlab figure
            %     FEVD figure

            % periods to plot
            fevd_periods = size(fevd,1);
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88]; 
            for period=1:fevd_periods
                % create custom bar plot
                left_x_edge = period - 0.4;
                right_x_edge = period + 0.4;
                median = fevd(period,1);
                lower = fevd(period,2);
                upper = fevd(period,3);
                hold on;
                x_patch = [left_x_edge, right_x_edge, right_x_edge, left_x_edge];
                y_patch = [upper upper lower lower];
                fevd_patch = fill(x_patch, y_patch, 'g');
                set(fevd_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
                plot([left_x_edge right_x_edge], [median median], 'LineWidth', 2, 'color', [0 0.5 0]);
                plot([left_x_edge left_x_edge], [lower upper], 'LineWidth', 0.5, 'color', [0 0 0]);
                plot([right_x_edge right_x_edge], [lower upper], 'LineWidth', 0.5, 'color', [0 0 0]);     
                plot([left_x_edge right_x_edge], [upper upper], 'LineWidth', 0.5, 'color', [0 0 0]);
                plot([left_x_edge right_x_edge], [lower lower], 'LineWidth', 0.5, 'color', [0 0 0]);
                hold off
            end
            % set graphic limits, background and font size for ticks
            hold on;
            plot([0 fevd_periods+1],[1 1], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([fevd_periods+1 fevd_periods+1],[0 1], 'LineWidth', 0.001, 'color', [0 0 0]);
            hold off
            set(gca,'XLim', [0 fevd_periods+1]);
            set(gca,'YLim', [0 1]);
            set(gca, 'color', [.9 .9 .9]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            % title
            title(['FEVD: ' name '\_' shock], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_fevd_joint(fevd, name, shocks, n)

            % var_fevd_joint(fevd, name, shocks, n)
            % produces FEVD figure for var model, single variable to all shocks
            % 
            % parameters:     
            % fevd: matrix of size (fevd_periods,n)
            %     fevd values for all shocks
            % name: char
            %     name of variable for which figure is produced
            % shocks: str array
            %     list of structural shocks for which figure is produced         
            % n: int
            %     number of endogenous variables
            %     
            % returns:
            % fig: matlab figure
            %     FEVD figure   

            % periods to plot
            fevd_periods = size(fevd,1);
            cum_fevd = cumsum(fevd,2);
            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            colors = gu.make_colors(n);
            for period=1:fevd_periods
                % create custom bar plot
                left_x_edge = period - 0.4;
                right_x_edge = period + 0.4;
                x_patch = [left_x_edge, right_x_edge, right_x_edge, left_x_edge];
                for i=1:n
                    if i == 1
                        lower = 0;
                    else
                        lower = cum_fevd(period,i-1);
                    end
                    upper = cum_fevd(period,i);
                    y_patch = [upper upper lower lower];
                    hold on;
                    if period == 1 && i < 16
                        legend('AutoUpdate','on');
                        fevd_patch = fill(x_patch, y_patch, 'g', 'DisplayName', char(shocks(i)));
                        legend('AutoUpdate','off');
                    else
                        fevd_patch = fill(x_patch, y_patch, 'g');
                    end
                    set(fevd_patch, 'FaceColor', colors(i,:), 'edgecolor', 'none', 'facealpha', 0.4);
                    plot([left_x_edge right_x_edge], [upper upper], 'LineWidth', 0.8, 'color', [0 0 0]);
                    hold off
                end
                hold on
                plot([left_x_edge left_x_edge], [0 1], 'LineWidth', 0.8, 'color', [0 0 0]);
                plot([right_x_edge right_x_edge], [0 1], 'LineWidth', 0.8, 'color', [0 0 0]);
                hold off
            end
            % set graphic limits, background and font size for ticks
            hold on;
            plot([0 fevd_periods+1],[1 1], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([fevd_periods+1 fevd_periods+1],[0 1], 'LineWidth', 0.001, 'color', [0 0 0]);
            hold off
            set(gca,'XLim', [0 fevd_periods+1]);
            set(gca,'YLim', [0 1]);
            set(gca, 'color', [.9 .9 .9]);
            set(gca, 'FontSize',13);
            % set figure and plot background color
            set(0,'defaultfigurecolor',[1 1 1]);
            grid on;
            set(gca, 'GridColor', [.3 .3 .3]);
            % create top and right axes
            box off;
            % title
            title(['FEVD: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            legend('Location','southoutside','Orientation','horizontal','NumColumns',5);
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_fevd_all(fevd, variables, shocks, n)

            % var_fevd_joint(fevd, name, shocks, n)
            % produces FEVD figure for var model, all variables to all shocks
            % 
            % parameters:     
            % fevd: matrix of size (n,n,fevd_periods,3)
            %     fevd values for all shocks
            % variables: str array
            %     name of variables for which figure is produced
            % shocks: char
            %     name of shocks for which figure is produced        
            % n: int
            %     number of endogenous variables
            %     
            % returns:
            % fig: matlab figure
            %     FEVD figure

            % periods to plot
            fevd_periods = size(fevd,1);
            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            x_patch = [fevd_periods;flipud(fevd_periods)];
            colors = gu.make_colors(n);
            % loop over variables
            for i=1:n
                axes('Units', 'normalized', 'Position', positions(i,:));
                cum_fevd = cumsum(fevd(:,:,i,1),2);
                for period=1:fevd_periods
                    % create custom bar plot
                    left_x_edge = period - 0.4;
                    right_x_edge = period + 0.4;
                    x_patch = [left_x_edge, right_x_edge, right_x_edge, left_x_edge];
                    for i=1:n
                        if i == 1
                            lower = 0;
                        else
                            lower = cum_fevd(period,i-1);
                        end
                        upper = cum_fevd(period,i);
                        y_patch = [upper upper lower lower];
                        hold on;
                        if period == 1 && i < 16
                            legend('AutoUpdate','on');
                            fevd_patch = fill(x_patch, y_patch, 'g', 'DisplayName', char(shocks(i)));
                            legend('AutoUpdate','off');
                        else
                            fevd_patch = fill(x_patch, y_patch, 'g');
                        end
                        set(fevd_patch, 'FaceColor', colors(i,:), 'edgecolor', 'none', 'facealpha', 0.4);
                        plot([left_x_edge right_x_edge], [upper upper], 'LineWidth', 0.8, 'color', [0 0 0]);
                        hold off
                    end
                    hold on
                    plot([left_x_edge left_x_edge], [0 1], 'LineWidth', 0.8, 'color', [0 0 0]);
                    plot([right_x_edge right_x_edge], [0 1], 'LineWidth', 0.8, 'color', [0 0 0]);
                    hold off
                end
                % set graphic limits, background and font size for ticks
                hold on;
                plot([0 fevd_periods+1],[1 1], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([fevd_periods+1 fevd_periods+1],[0 1], 'LineWidth', 0.001, 'color', [0 0 0]);
                hold off
                set(gca,'XLim', [0 fevd_periods+1]);
                set(gca,'YLim', [0 1]);
                set(gca, 'color', [.9 .9 .9]);
                set(gca, 'FontSize',13);
                % set figure and plot background color
                set(0,'defaultfigurecolor',[1 1 1]);
                grid on;
                set(gca, 'GridColor', [.3 .3 .3]);
                % create top and right axes
                box off;
                % title
                name = char(variables(i));
                title(['FEVD: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                legend('Location','southoutside','Orientation','horizontal','NumColumns',5);
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [fig] = var_hd_single_variable(hd, name, shock, dates)

            % var_hd_single_variable(hd, name, shock, dates)
            % produces HD figure for var model, single variable
            % 
            % parameters:     
            % hd: numpy ndarray of size (T,3)
            %     hd values, median, lower and upper bounds
            % name: char
            %     name of variable for which figure is produced
            % shock: char
            %     name of shock for which figure is produced  
            % dates: datetime array of size (T)
            %     index of in-sample dates            
            %     
            % returns:
            % fig: matlab figure
            %     HD figure

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            [min_YLim max_YLim] = gu.set_min_and_max(hd, 0.07, 0.1);
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            % plot irf with patched credibility intervals
            hold on;
            x_patch = [dates;flipud(dates)];
            y_patch = [hd(:,2);flipud(hd(:,3))];
            fitted_patch = fill(x_patch, y_patch, 'g');
            set(fitted_patch, 'FaceColor', [0.3 0.7 0.2], 'edgecolor', 'none', 'facealpha', 0.4);
            plot([dates(1) dates(end)], [0 0], 'LineStyle', '--', 'color', [0 0 0], 'LineWidth', 1);
            plot(dates, hd(:,1), 'LineWidth', 1.3, 'color', [0 0.5 0]);
            plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
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
            % title
            title(['HD: ' name '\_' shock], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_hd_joint(hd, name, shocks, dates, n, T)

            % var_hd_joint(hd, name, shocks, dates, n, T)
            % produces HD figure for var model, single variable
            % 
            % parameters:     
            % hd: matrix of size (T,n)
            %     hd valuesfor all shocks
            % name: char
            %     name of variable for which figure is produced
            % shocks: str array
            %     name of shocks for which figure is produced         
            % dates: datetime index of size (T)
            %     index of in-sample dates
            % n: int
            %     number of endogenous variables
            % T: int
            %     number of sample periods
            %     
            % returns:
            % fig: matlab figure
            %     HD figure   

            % create figure
            fig = figure('Position', [100 100 660 470], 'Visible', 'off');
            ax = gca;
            ax.Position = [0.08 0.06 0.88 0.88];
            x_patch = [dates;flipud(dates)];
            colors = gu.make_colors(n);
            % trend
            cum_hd = sum(hd,2);
            hold on
            legend('AutoUpdate','on');
            plot(dates, cum_hd, 'LineWidth', 1.5, 'color', [0 0 0], 'DisplayName', 'trend');
            legend('AutoUpdate','off');
            hold off
            % positive contributions
            positive_hd = hd;
            positive_hd(positive_hd<0) = 0;
            cum_positive_hd = cumsum(positive_hd,2);
            cum_positive_hd = [zeros(T,1) cum_positive_hd]; 
            for i=1:n
                y_patch = [cum_positive_hd(:,i+1); flipud(cum_positive_hd(:,i))];
                hold on;
                if i < 15
                    legend('AutoUpdate','on');
                    fitted_patch = fill(x_patch, y_patch, 'g', 'DisplayName', char(shocks(i)));
                    legend('AutoUpdate','off');
                else
                    fitted_patch = fill(x_patch, y_patch, 'g');
                end
                set(fitted_patch, 'FaceColor', colors(i,:), 'edgecolor', 'none', 'facealpha', 0.4);
                plot(dates, cum_positive_hd(:,i+1), 'LineWidth', 0.01, 'color', [0.6 0.6 0.6]);
                hold off
            end
            % negative contributions
            negative_hd = hd;
            negative_hd(positive_hd>0) = 0;
            cum_negative_hd = cumsum(negative_hd,2);
            cum_negative_hd = [zeros(T,1) cum_negative_hd]; 
            for i=1:n
                y_patch = [cum_negative_hd(:,i+1); flipud(cum_negative_hd(:,i))];
                hold on
                fitted_patch = fill(x_patch, y_patch, 'g');
                set(fitted_patch, 'FaceColor', colors(i,:), 'edgecolor', 'none', 'facealpha', 0.4);
                plot(dates, cum_negative_hd(:,i+1), 'LineWidth', 0.01, 'color', [0.6 0.6 0.6]);
                hold off
            end
            % set graphic limits, background and font size for ticks
            [min_YLim max_YLim] = gu.set_min_and_max([cum_positive_hd cum_negative_hd], 0.07, 0.1);
            hold on
            plot(dates, cum_hd, 'LineWidth', 1.5, 'color', [0 0 0]);
            plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
            hold off
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
            % title
            title(['HD: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
            legend('Location','southoutside','Orientation','horizontal','NumColumns',5);
            % command to preserve grey background color when saving as image
            fig.InvertHardcopy = 'off';
        end


        function [fig] = var_hd_all(hd, variables, shocks, dates, n, T)

            % var_hd_all(hd, variables, shocks, dates, n, T)
            % produces HD figure for var model, all variables to all shocks
            % 
            % parameters:     
            % hd: matrix of size (n,n,T,3)
            %     hd values for all shocks
            % variables: str array
            %     name of variables for which figure is produced
            % shocks: str array
            %     name of shocks for which figure is produced         
            % dates: datetime index of size (T)
            %     index of in-sample dates
            % n: int
            %     number of endogenous variables
            % T: int
            %     number of sample periods
            %     
            % returns:
            % fig: matlab figure
            %     HD figure 

            % get plot dimensions
            columns = ceil(n ^ 0.5);
            rows = columns;
            % create figure
            fig = figure('Position', [100 100 660*columns 470*rows], 'Visible', 'off');
            [positions] = gu.make_subplot_positions(n, rows, columns, 0.07, 0.05, 0.04, 0.05);
            x_patch = [dates;flipud(dates)];
            colors = gu.make_colors(n);
            % loop over variables
            for i=1:n
                % initiate subplot
                axes('Units', 'normalized', 'Position', positions(i,:));
                % recover HD for this variable
                variable_hd = hd(:,:,i,1);
                % trend
                cum_hd = sum(variable_hd,2);
                hold on
                legend('AutoUpdate','on');
                plot(dates, cum_hd, 'LineWidth', 1.5, 'color', [0 0 0], 'DisplayName', 'trend');
                legend('AutoUpdate','off');
                hold off
                % positive contributions
                positive_hd = variable_hd;
                positive_hd(positive_hd<0) = 0;
                cum_positive_hd = cumsum(positive_hd,2);
                cum_positive_hd = [zeros(T,1) cum_positive_hd]; 
                for j=1:n
                    y_patch = [cum_positive_hd(:,j+1); flipud(cum_positive_hd(:,j))];
                    hold on;
                    if j < 15
                        legend('AutoUpdate','on');
                        fitted_patch = fill(x_patch, y_patch, 'g', 'DisplayName', char(shocks(j)));
                        legend('AutoUpdate','off');
                    else
                        fitted_patch = fill(x_patch, y_patch, 'g');
                    end
                    set(fitted_patch, 'FaceColor', colors(j,:), 'edgecolor', 'none', 'facealpha', 0.4);
                    plot(dates, cum_positive_hd(:,j+1), 'LineWidth', 0.01, 'color', [0.6 0.6 0.6]);
                    hold off
                end
                % negative contributions
                negative_hd = variable_hd;
                negative_hd(positive_hd>0) = 0;
                cum_negative_hd = cumsum(negative_hd,2);
                cum_negative_hd = [zeros(T,1) cum_negative_hd]; 
                for j=1:n
                    y_patch = [cum_negative_hd(:,j+1); flipud(cum_negative_hd(:,j))];
                    hold on
                    fitted_patch = fill(x_patch, y_patch, 'g');
                    set(fitted_patch, 'FaceColor', colors(j,:), 'edgecolor', 'none', 'facealpha', 0.4);
                    plot(dates, cum_negative_hd(:,j+1), 'LineWidth', 0.01, 'color', [0.6 0.6 0.6]);
                    hold off
                end
                % set graphic limits, background and font size for ticks
                [min_YLim max_YLim] = gu.set_min_and_max([cum_positive_hd cum_negative_hd], 0.07, 0.1);
                hold on
                plot(dates, cum_hd, 'LineWidth', 1.5, 'color', [0 0 0]);
                plot([dates(1) dates(end)],[max_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                plot([dates(end) dates(end)],[min_YLim max_YLim], 'LineWidth', 0.001, 'color', [0 0 0]);
                hold off
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
                % title
                name = char(variables(i));
                title(['HD: ' name], 'FontWeight', 'bold','FontSize',14,'FontName','Serif');
                legend('Location','southoutside','Orientation','horizontal','NumColumns',5);
                % command to preserve grey background color when saving as image
                fig.InvertHardcopy = 'off';
            end
        end


        function [min_YLim max_YLim] = set_min_and_max(data, min_space, max_space)
            
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
  

        function show_and_save(current_figure, show, save, path, file_name)
            
            % show_and_save(current_figure, show, save, path, file_name)
            % display and save a given figure, if requested
            % 
            % parameters:
            % current_figure : matlab figure
            %     figure to display and save
            % show : bool
            %     if true, display the figure
            % save : bool
            %     if true, save figure as image
            % path : str
            %     path to folder where figure is saved
            % file_name : str
            %     name of file to save the figure
            %         
            % returns:
            % none
            
            % save in folder as image if requested
            if save
                full_path = fullfile(path, file_name);
                saveas(current_figure, full_path);
            end
            % display figure if requested
            if show
                figure(current_figure);
            end
        end


        function [colors] = make_colors(n)
            
            % make_colors(n)
            % return an array of n color triplets, in a determined order
            % 
            % parameters:
            % n : int
            %     number of color triplets to return
            %         
            % returns:
            % colors: ndarray of size(n,3)
            %     array of column triplets 
            
            % define fixed colors
            fixed_colors = [
                0.196 0.804 0.196;  % 'limegreen'
                1.000 0.647 0.000;  % 'orange'
                0.000 1.000 1.000;  % 'cyan'
                1.000 0.000 0.000;  % 'red'
                0.580 0.000 0.827;  % 'darkviolet'
                0.647 0.165 0.165;  % 'brown'        
                1.000 0.714 0.757;  % 'lightpink' 
                0.753 0.753 0.753;  % 'silver'    
                0.502 0.502 0.000;  % 'olive'    
                0.000 0.000 1.000;  % 'blue'   
                0.804 0.521 0.247;  % 'peru'   
                0.541 0.169 0.886;  % 'blueviolet'   
                1.000 1.000 0.000;  % 'yellow' 
                0.498 1.000 0.000;  % 'chartreuse' 
                0.863 0.078 0.235;  % 'crimson' 
                0.688 0.389 0.135;  %  random colors from here on
                0.721 0.525 0.310;
                0.486 0.889 0.934;
                0.358 0.572 0.322;
                0.594 0.338 0.392;
                0.890 0.227 0.623;
                0.084 0.833 0.787;
                0.239 0.876 0.059;
                0.336 0.150 0.450;
                0.796 0.231 0.052;
                0.405 0.199 0.091;
                0.580 0.299 0.672;
                0.200 0.942 0.365;
                0.105 0.629 0.927;
                0.440 0.955 0.500;
                0.425 0.620 0.995;
                0.949 0.460 0.758;
                0.497 0.529 0.786;
                0.415 0.734 0.711;
                0.932 0.115 0.729;
                0.927 0.968 0.015;
                0.864 0.981 0.957;
                0.149 0.973 0.890;
                0.822 0.480 0.232;
                0.802 0.924 0.266;
                0.539 0.443 0.931;
                0.041 0.732 0.614;
                0.028 0.719 0.016;
                0.758 0.513 0.929;
                0.066 0.841 0.067;
                0.344 0.430 0.966;
                0.562 0.259 0.242;
                0.888 0.226 0.125;
                0.288 0.586 0.554;
                0.810 0.560 0.288;
                0.413 0.818 0.627;
                0.959 0.369 0.553;
                0.594 0.848 0.145;
                0.407 0.910 0.043;
                0.823 0.415 0.830;
                0.010 0.365 0.079;
                0.653 0.274 0.703;
                0.944 0.127 0.865;
                0.059 0.381 0.430;
                0.489 0.976 0.776;
                0.309 0.270 0.863;
                0.881 0.511 0.344;
                0.995 0.316 0.183;
                0.880 0.812 0.668;
                0.958 0.926 0.748;
                0.861 0.247 0.141;
                0.670 0.715 0.167;
                0.396 0.910 0.561;
                0.578 0.194 0.526;
                0.523 0.089 0.982;
                0.571 0.006 0.773;
                0.978 0.590 0.320;
                0.188 0.673 0.195;
                0.578 0.602 0.962;
                0.072 0.500 0.744;
                0.177 0.388 0.063;
                0.726 0.088 0.395;
                0.874 0.472 0.913;
                0.766 0.915 0.127;
                0.074 0.070 0.869;
                0.634 0.497 0.164;
                0.674 0.318 0.711;
                0.460 0.507 0.790;
                0.093 0.579 0.197;
                0.808 0.489 0.989;
                0.183 0.963 0.801;
                0.481 0.814 0.603;
                0.655 0.914 0.065;
                0.835 0.382 0.326;
                0.994 0.781 0.486;
                0.423 0.878 0.087;
                0.708 0.789 0.799;
                0.322 0.797 0.225;
                0.362 0.417 0.541;
                0.113 0.407 0.000;
                0.744 0.852 0.139;
                0.704 0.821 0.982;
                0.844 0.424 0.980;
                0.974 0.504 0.753
                ]; 
            if n > 100
                rng(0);
                random_colors = reshape(rand(3*(n-100),1),[3 n-100])';
                rng("shuffle");
                colors = [fixed_colors; random_colors];
            else
                colors = fixed_colors(1:n,:);
            end
        end


        function [positions] = make_subplot_positions(subplots, rows, columns, ...
                 row_space, column_space, row_margin, column_margin)

            % make_subplot_positions(subplots, rows, columns, ...
            %     ... row_space, column_space, row_margin, column_margin)
            % compute positions of subplots with tight spaces
            % 
            % parameters:     
            % subplots: int
            %     total number of graphs to plot as subplots
            % rows: int
            %     number of subplot rows in the subplot grid
            % columns: int
            %     number of subplot columns in the subplot grid
            % row_space: float between 0 and 1
            %     row (vertical) space between two subplots, as proportion of figure
            % column_space: float between 0 and 1
            %     column (horizontal) space between two subplots, as proportion of figure
            % row_margin: float between 0 and 1
            %     row (vertical) space on figure margins, as proportion of figure
            % column_margin: float between 0 and 1
            %     column (horizontal) space on figure margins, as proportion of figure            
            %     
            % returns:
            % fig: matlab figure
            %     fitted figure

            subplot_height = (1 - 2 * row_margin - (rows - 1) * row_space) / rows; 
            subplot_width = (1 - column_margin - 0.02 - (columns - 1) * column_space) / columns; 
            y_position = 1 - row_margin - subplot_height;
            i = 0;
            for j = 1:rows
                x_position = column_margin;
                for k = 1:columns
                    i = i + 1;
                    if i <= subplots
                        positions(i,:) = [x_position y_position subplot_width subplot_height];
                        x_position = x_position + subplot_width + column_space;
                    end
                end
            y_position = y_position - subplot_height - row_space;
            end
        end

    end
        
end
         