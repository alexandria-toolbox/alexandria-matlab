classdef cu
    

    % cu stands for console utilities
    % A class containing static methods for console display

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function progress_bar(iteration, total_iterations, application_tag)
            
            % progress_bar(iteration, total_iterations, application_tag)
            % progress bar for MCMC algorithms
            %
            % parameters:
            % iteration: int
            %   current iteration of MCMC algorithm
            % total_iterations: int
            %   total number of iterations of MCMC algorithm 
            % application_tag: char
            %   string flagging the application being estimated by MCMC
            %
            % returns:
            % none

            if iteration == 1
                % print application tag (which application is being run)
                fprintf('%s\n', application_tag);
                % print iteration 1 progress bar
                string = [repmat(' ', 1, numel(num2str(total_iterations))) ...
                     '1/' num2str(total_iterations) '  [' repmat('.', 1, 33) ...
                     ']   0%'];
                 fprintf(string);
            elseif iteration == total_iterations
                % length of string to overwrite with new string
                string_length = 43 + 2 * numel(num2str(total_iterations));
                % print final iteration progress bar
                string = [repmat('\b', 1, string_length) ...
                    num2str(total_iterations) '/' num2str(total_iterations) ...
                    '  [' repmat('=',1,33) ']  —  done'];
                fprintf(string);
                fprintf('\n');
            else
                % integer defining the arrow position in the progress bar
                arrow_position = floor(iteration * 33.3333 / total_iterations);
                % integer percentage of progress
                percentage = floor(iteration * 100 / total_iterations);
                string_percentage = [repmat(' ', 1, double(percentage < 10)) ...
                    num2str(percentage) '%'];
                % length of string to overwrite with new string
                string_length = 43 + 2 * numel(num2str(total_iterations));
                string = [repmat('\b', 1, string_length) ...
                    repmat(' ', 1, numel(num2str(total_iterations)) ...
                    - numel(num2str(iteration))) num2str(iteration) '/' ...
                    num2str(total_iterations) '  [' ...
                    repmat('=', 1, arrow_position - 1) ...
                    repmat('>', 1, double(arrow_position ~= 0)) ...
                    repmat('.', 1, 33 - arrow_position) ']  ' ...
                    string_percentage '%'];
                fprintf(string);
                % pause(1e-08);
            end
        end
        
        
        function progress_bar_complete(application_tag)
    
            % progress_bar_complete(application_tag)
            % pseudo progress bar used to show completion when no actual MCMC is run
            %
            % parameters:
            % application_tag: char
            %   string flagging the application being estimated by MCMC
            %
            % returns:
            % none

            % print application tag (which application is being run)
            fprintf('%s\n',application_tag);   
            % print pseudo progress bar
            string = ['  — / —    [' repmat('=',1,33) ']  —  done'];
            fprintf('%s\n', string);
        end
        
    
        function optimization_completion(success)
        
            % optimization_completion(success)
            % binary message that indicates success or failure of optimization
            % 
            % parameters:      
            % success: bool
            %     boolean indicating whether optimization has succeeded
            %     
            % returns:
            % none

            if success
                fprintf('Optimization conducted succesfully.\n');
            else
                fprintf('Warning! Optimization failed. Prior values may be unreliable.\n')
            end
        end
        
        
        function print_message(message)
            
            % print_message(message)
            % a simple print function that prints any message
            %
            % parameters:       
            % message: char
            %     the message to print
            %
            % returns:
            % none   

            disp(message);
        end
        
        
        function print_message_to_complete(message)
            
            % print_message_to_complete(message)
            % a print function that prints a message that will be completed by the next message
            %
            % parameters:       
            % message: char
            %     the message to print
            %
            % returns:
            % none   

            fprintf(message);
        end        
        
        
        function print_string_array(string_array)
        
            % print_string_list(string_list)
            % a function that prints an array of strings, line by line
            % 
            % parameters:        
            % string_array: str array
            %     array containing the strings to print
            % 
            % returns
            % none
            
            for i=1:numel(string_array)
                disp(string_array(i));
            end
        end
        

        function write_string_array(string_array, filepath)

            % print_string_list(string_list)
            % a function that writes a list of strings, line by line, in a specified path
            % 
            % parameters      
            % string_list: str array
            %     string array containing the strings to print
            % filepath: char
            %     full path to file (must include file name)
            % 
            % returns
            % none
            
            file = fopen(filepath, 'w');
            for i=1:numel(string_array)
                fprintf(file, '%s\n', string_array(i));
            end
            fclose(file);
        end
            

        function check_path(path)
        
            % check_path(path)
            % checks whether folder given by path exists, and create it if needed
            % 
            % parameters:     
            % path: char
            %     path to folder
            %     
            % returns:
            % none 
        
            if ~(exist(path) == 7)
                mkdir(path)
            end
        end

 
        function [header] = alexandria_header()

            % [header] = alexandria_header()
            % a function that creates a list of strings constituting Alexandria's header
            % 
            % parameters      
            % none
            % 
            % returns
            % header: str array
            %     the list of strings constituting Alexandria's header
     
            header = [...
                      "                                                                               "; ...
                      "      ======================================================================== "; ...
                      "         _____   __                                   _                        "; ...
                      "        / __  / / / ___  __  __   ____    ____    ___/ / ____  (_) ____        "; ...
                      "       / /_/ / / / / _ \ \ \/ /  / __ `/ / __ \ / __  / / __/ / / / __ `/      "; ...
                      "      / __  / / / /  __/  |  |  / /_/ / / / / // /_/ / / /   / / / /_/ /       "; ...
                      "     /_/ /_/ /_/  \___/  /_/\_\ \__,_/ /_/ /_/ \____/ /_/   /_/  \__,_/        "; ...
                      "                                                                               "; ...
                      "     The library of Bayesian time-series models                                "; ...
                      "     V 1.0 - Copyright Ⓒ  Romain Legrand                                      "; ...
                      "   ========================================================================    "; ...
                      "                                                                               "; ...
                      "                                                                               "; ...
                      ];
        end
        

        function print_alexandria_header()
            
            % print_alexandria_header()
            % display Alexandria header on console
            % 
            % parameters:     
            % none
            %     
            % returns:
            % none      
            
            % get Alexandria header
            header = cu.alexandria_header();
            % print header
            cu.print_string_array(header);
        end


        function print_start_message()
        
            % print_start_message()
            % display start message
            % 
            % parameters:     
            % none
            %     
            % returns:
            % none  
            
            cu.print_message('Starting estimation of your model...');
            cu.print_message(' ');
        end


        function print_completion_message(progress_bar)
        
            % print_completion_message(progress_bar)
            % display completion message
            % 
            % parameters:     
            % progress_bar: bool
            %     if yes, progress bar is displayed
            %     
            % returns:
            % none  
            
            if progress_bar
                cu.print_message(' ');
            end
            cu.print_message('Estimation completed successfully.');
            cu.print_message(' ');
            cu.print_message(' ');
        end


        function [formatted_number] = format_number(number)
            
            % [formatted_number] = format_number(number)
            % formats any number into a 10-character string
            %
            % parameters       
            % number: float
            %     number to be formatted
            %
            % returns
            % formatted_number: char
            %     string containing the formatted number

            % if number is exactly 0, format it normally
            if number == 0
                formatted_number = '     0.000';    
            % if number is of regular length, use decimal notation with 3 decimals
            elseif 0.0001 < abs(number) && abs(number) < 100000
                formatted_number = sprintf('%10.3f', number);
            else
            % if value is too small or too large, switch to exponential notation    
                formatted_number = sprintf('%10.3e', number);
            end
        end
        
        
        function [formatted_name] = format_name(name)
            
            % [formatted_name] = format_name(name)
            % formats any variable name as a string containing a maximum of 20 characters
            %
            % parameters    
            % name: char
            %     variable name to be formatted
            %
            % returns
            % formatted_name: char
            %     string containing the formatted name

            % if name is less than 20 characters, return it left-justified
            formatted_name = sprintf('%-20s', cu.shorten_string(name, 20));
        end
        
        
        function [string] = shorten_string(string, n)
            
            % [string] = shorten_string(string, n)
            % allows for a maximum of n characters in string; if longer, string is shortened with '..'
            %
            % parameters       
            % string: char
            %     string to be formatted
            %
            % returns
            % string: char
            %     string containing at most n characters

            if numel(string) > n
                string = [string(1:n-2) '..'];
            end
        end
        
        
        function [line] = equal_dashed_line()

            % [line] = equal_dashed_line()
            % return a line of 80 equal signs ('=')
            %
            % parameters     
            % none
            %
            % returns
            % line: char
            %     string containing the line
            
            line = convertCharsToStrings(repmat('=', 1, 80));
        end
        
        
        function [line] = hyphen_dashed_line()

            % [line] = hyphen_dashed_line()
            % return a line of 80 hyphen signs ('-')
            %
            % parameters     
            % none
            %
            % returns
            % line: char
            %     string containing the line
            
            line = convertCharsToStrings(repmat('-', 1, 80));
        end        
        
        
        function [header] = model_header(model)

            % [header] = model_header(model)
            % return a results header with model name and two wrapping equal dashed lines
            %
            % parameters     
            % model: char
            %     char containing the model name
            %
            % returns
            % header: str array
            %     array of string containing the model header

            % initiate header, add equal dashed line
            header = [cu.equal_dashed_line()];
            % center-justify model name
            header = [header; cu.center_align(model, 80)];
            % add final equal dashed line
            header = [header; cu.equal_dashed_line()];  
        end
        
        
        function [centered_string] = center_align(string, n)
            
            % [centered_string] = center_align(string)
            % returns a n-character char with the string center-aligned
            %
            % parameters:
            % string: char
            %     the char to be center-aligned
            %
            % returns
            % centered_string: char
            %     n-character char with the string center-aligned
            
            string_length = numel(string);
            left_pad = repmat(' ', 1, floor((n - string_length) / 2));
            right_pad = repmat(' ', 1, ceil((n - string_length) / 2));
            centered_string = [left_pad string right_pad];
        end
           
        
        function [header] = equation_header(string)
            
            % equation_header(string)
            % return a header with string, and two wrapping dashed lines
            % 
            % parameters:      
            % string: char
            %     string containing the model name
            %     
            % returns:
            % header: str array
            %     list of string containing the model header
            
            % initiate header, add equal dashed line
            header = [cu.equal_dashed_line()];
            % center-justify model name
            header = [header; cu.center_align(string, 80)];
            % add final equal dashed line
            header = [header; cu.hyphen_dashed_line()];  
        end


        function [header] = intermediate_header(string)
            
            % intermediate_header(string)
            % return an intermediate header with string, and two wrapping dashed lines
            % 
            % parameters:      
            % string: char
            %     string containing the model name
            %     
            % returns:
            % header: str array
            %     list of string containing the model header
            
            % initiate header, add equal dashed line
            header = [cu.hyphen_dashed_line()];
            % center-justify model name
            header = [header; cu.center_align(string, 80)];
            % add final equal dashed line
            header = [header; cu.hyphen_dashed_line()];  
        end


        function [regressors] = make_regressors(endogenous, exogenous, constant, trend, quadratic_trend, n, p)

            % make_regressors(endogenous, exogenous, constant, trend, quadratic_trend, n, p)
            % return a string array of regressors
            % 
            % parameters:      
            % endogenous: str array
            %     list of endogenous variables
            % exogenous: str array
            %     list of exogenous variables   
            % constant: bool
            %     if true, a constant is added to the model
            % trend: bool
            %     if true, a trend is added to the model
            % quadratic_trend: bool
            %     if true, a quadratic trend is added to the model
            % n : int
            %     number of endogenous variables
            % p : int
            %     number of lags
            %     
            % returns:
            % regressors: str array
            %     list of string containing the model regressors

            regressors = string([]);
            if constant
                regressors = [regressors 'constant'];
            end
            if trend
                regressors = [regressors 'trend'];
            end
            if quadratic_trend
                regressors = [regressors 'quadratic trend'];
            end
            if ~isequal(exogenous,"none")
                regressors = [regressors exogenous];
            end
            for i=1:n
                for j=1:p
                    regressors = [regressors [char(endogenous(i)) ' (-' num2str(j) ')']];
                end
            end
        end


        function [index] = make_index(n, m, p, k)

            % make_index(n, m, p, k)
            % return an array of indices for VAR coefficients
            % 
            % parameters:      
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables            
            % p : int
            %     number of lags
            % k : int
            %     number of coefficients per VAR equation
            %     
            % returns:
            % index: matrix of size (k,1)
            %     array of indices

            index = zeros(k,1);
            i = 0;
            for j=1:m
                i = i + 1;
                index(i) = j;
            end
            for g=1:n
                for h=1:p
                    i = i + 1;
                    index(i) = m + (h-1) * n + g;
                end
            end
        end


        function [line] = variance_line(residual_variance, shock_variance)
            
            % variance_line(residual_variance, shock_variance)
            % return a line with residual and shock variance estimate
            % 
            % parameters:      
            % variable: str
            %     string containing the model name
            %     
            % returns:
            % residual_variance: float
            %     residual variance estimate
            % residual_variance: float or empty char
            %     shock variance estimate    
            
            formatted_residual_variance = cu.format_number(residual_variance);
            left_element = ['residual variance:' sprintf('%+20s', formatted_residual_variance)];
            if isnumeric(shock_variance)
                formatted_shock_variance = cu.format_number(shock_variance);
                right_element = ['shock variance:' sprintf('%+23s', formatted_shock_variance)];
            else
                right_element = repmat(' ', 1, 38);
            end
            line = [left_element '    ' right_element];
        end


        function [lines] = variance_covariance_summary(Sigma, n, endogenous_variables, tag)
        
            % variance_covariance_summary(Sigma, n, endogenous_variables)
            % return a set of lines that summarizes a variance-covariance matrix
            % 
            % parameters:      
            % Sigma: ndarray of shape (n,n)
            %     spd variance-covariance matrix
            % n : int
            %     number of endogenous variables
            % endogenous_variables: str array
            %     list of endogenous variables
            % tag: str
            %     string providing the command to use to display full matrix
            %     
            % returns:
            % lines: str list
            %     string containing variance-covariance matrix summary   
        

            lines = [];
            dimension = min(size(Sigma,1), 6);
            header_line = repmat(' ', 1, 10);
            for i=1:dimension
                header_line = [header_line sprintf('%+11s', cu.shorten_string(char(endogenous_variables(i)), 10))];
            end
            header_line = cu.string_line(header_line);
            lines = [lines;header_line];
            for i=1:dimension
                current_line = sprintf('%-10s', cu.shorten_string(char(endogenous_variables(i)), 10));
                for j=1:dimension
                    current_line = [current_line ' ' cu.format_number(Sigma(i,j))];
                end
                if n>6
                    current_line = [current_line ' ...'];
                end
                current_line = cu.string_line(current_line);
                lines = [lines;current_line];
            end
            if n>6
                current_line = '  ⋮               ⋮          ⋮          ⋮          ⋮          ⋮          ⋮    ⋱ ';
                lines = [lines;current_line];
                lines = [lines;repmat(' ', 1, 80)];
                current_line = ['output is too long, use ' tag ' to obtain full view'];
                lines = [lines;cu.string_line(current_line)];
            end
        end


        function [line] = forecast_evaluation_line(variable, value_1, value_2, value_3, value_4, value_5)
            
            % forecast_evaluation_line(variable, value_1, value_2, value_3, value_4, value_5)
            % return a line with variable name and formatted numerical values
            % 
            % parameters:  
            % variable: char
            %     variable name
            % value_1: float
            %     first value to display on line
            % value_2: float
            %     second value to display on line
            % value_3: float
            %     third value to display on line
            % value_4: float
            %     fourth value to display on line
            % value_5: float
            %     fifth value to display on line
            %   
            % returns:
            % line: str
            %     string containing formatted forecast evaluation criteria summary         
            
            line = sprintf('%-10s', cu.shorten_string(variable, 10));
            line = [line ' ' cu.format_number(value_1)];
            line = [line ' ' cu.format_number(value_2)];
            line = [line ' ' cu.format_number(value_3)];
            line = [line ' ' cu.format_number(value_4)];
            line = [line ' ' cu.format_number(value_5)];
            line = cu.string_line(line);
        end


        function [lines] = forecast_evaluation_summary(log_score, joint_log_score, endogenous_variables, tag)
            
            % forecast_evaluation_summary(log_score, joint_log_score, endogenous_variables, tag)
            % return a set of lines that summarizes Bayesian forecast evaluation criteria
            %
            % parameters:      
            % log_score: matrix of size (forecast_periods,n)
            %     array of log score values
            % joint_log_score: matrix of size (forecast_periods,)
            %     array of joint log score values
            % endogenous_variables: str array
            %     list of endogenous variables
            % tag: char
            %     string providing the command to use to display full matrix
            %     
            % returns:
            % line: char
            %     list of strings containing formatted forecast evaluation criteria summary    
            
            lines = [];
            periods = size(log_score,1);
            dimension = min(periods, 6);
            header_line = repmat(' ', 1, 10);
            for i=1:dimension
                header_line = [header_line sprintf('%+11s', ['(+' num2str(i) ')'])];
            end
            if dimension == 6
                header_line = [header_line ' ...'];
            else
                header_line = [header_line '      (all)'];
            end
            header_line = cu.string_line(header_line);
            lines = [lines;header_line];
            for i=1:numel(endogenous_variables)
                line = sprintf('%-10s', cu.shorten_string(char(endogenous_variables(i)), 10));
                for j=1:dimension
                    value = log_score(j,i);
                    line = [line ' ' cu.format_number(value)];
                end
                if dimension == 6
                    line = [line ' ...'];
                else
                    value = joint_log_score(i);
                    line = [line ' ' cu.format_number(value)];
                end
                line = cu.string_line(line);
                lines = [lines;line];
            end
            if dimension == 6
                lines = [lines;repmat(' ', 1, 80)];
                current_line = ['output is too long, use ' tag ' to obtain full view'];
                lines = [lines;cu.string_line(current_line)];
            end
        end


        function [header] = coefficient_header(credibility_level)

            % [header] = coefficient_header(credibility_level)
            % return an estimation line with headers for coefficient, standard devations and credibility bounds
            %
            % parameters
            % credibility_level: float
            %     credibility level for model
            %
            % returns
            % line: char
            %     string containing coefficient header

            % calculate lower and upper bounds
            lower_bound = sprintf('%5.3f', 0.5 - credibility_level / 2);
            upper_bound = sprintf('%5.3f', 0.5 + credibility_level / 2);
            header = [repmat(' ', 1, 29) 'median' repmat(' ', 1, 8) ...
                      'std dev' repmat(' ', 1, 9) '[' lower_bound ...
                      repmat(' ', 1, 9) upper_bound ']'];
            % generate second line of header
            header = [header; cu.hyphen_dashed_line()];
        end
            
        
        function [line] = parameter_estimate_line(name, coefficient, standard_deviation, lower_bound, upper_bound)

            % [line] = parameter_estimate_line(name, coefficient, standard_deviation, lower_bound, upper_bound)
            % return an estimation line with formatted values for coefficient, standard devations and credibility bounds
            %
            % parameters
            % name: char
            %     name of variable to which parameter is related
            % coefficient: float
            %     coefficient value
            % standard_deviation: float (positive)
            %     standard deviation of parameter
            % lower_bound: float
            %     lower bound of parameter
            % upper_bound: float
            %     upper bound of parameter
            %
            % returns
            % line: char
            %     string containing coefficient summary

            line = [cu.format_name(name) repmat(' ', 1, 5) cu.format_number(coefficient) ...
                    repmat(' ', 1, 5) cu.format_number(standard_deviation) repmat(' ', 1, 5) ...
                    cu.format_number(lower_bound) repmat(' ', 1, 5) cu.format_number(upper_bound)];
        end
        
        
        function [line] = string_line(string)
            
            % string_line(name)
            % returns a line where name is left-aligned, and right padding is added to reach 80 characters
            %
            % parameters
            % string: char
            %     string to left-align
            %
            % returns
            % line: str
            %     80-character line with left-aligned string

            % left-justify the name, with pad to reach 80 characters
            line = sprintf('%-80s', string);
        end
     
            
        function [lines] = insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic)

            % [lines] = insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic)
            % returns the set of lines with the results for in-sample evaluation criteria
            %
            % parameters
            % ssr: float
            %     sum of squared residuals
            % r2: float
            %     coefficient of determination
            % adj_r2: float
            %     adjusted coefficient of determination
            % m_y: float
            %     log10 marginal likelihood
            % aic: float
            %     Akaike information criterion
            % bic: float
            %     Bayesian information criterion
            %
            % returns
            % lines: str array
            %     set of lines reporting in-sample evaluation criteria

            % initiate lines
            lines = [];
            % first row: r2 and ssr
            left_element = sprintf(['R2:' sprintf('%+35s', cu.format_number(r2))]);
            right_element = sprintf(['ssr:' sprintf('%+34s', cu.format_number(ssr))]);
            lines = [left_element  '    ' right_element];
            % second row: adjusted r2 and marginal likelihood
            left_element = sprintf(['adj. R2:' sprintf('%+30s', cu.format_number(adj_r2))]);
            if m_y
                right_element = sprintf(['log10 marg. lik.:' sprintf('%+21s', cu.format_number(m_y))]);
            else
                right_element = repmat(' ', 1, 38);
            end
            lines = [lines; [left_element  '    ' right_element]];
            % third row AIC and BIC
            if aic
                left_element = sprintf(['AIC:' sprintf('%+34s', cu.format_number(aic))]);
                right_element = sprintf(['BIC:' sprintf('%+34s', cu.format_number(bic))]);
                lines = [lines; [left_element  '    ' right_element]];
            end
        end
            

        function [lines] = forecast_evaluation_lines(rmse, mae, mape, theil_u, bias, log_score, crps)

            % [lines] = forecast_evaluation_lines(rmse, mae, mape, theil_u, bias, log_score, crps)
            % returns the set of lines with the results for out-of-sample evaluation criteria
            %
            % parameters
            % rmse: float
            %     root mean squared error for predictions
            % mae: float
            %     mean absolute error for predictions
            % mape: float
            %     mean absolute percentage error for predictions
            % theil_u: float
            %     Theil's U coefficient for predictions
            % bias: float
            %     bias for predictions
            % log_score: float
            %     log score for predictions
            % crps: float
            %     continuous rank probability score for predictions
            %
            % returns
            % lines: str array
            %     set of lines reporting out-of-sample evaluation criteria

            % first row: rmse
            left_element = sprintf(['rmse:' sprintf('%+33s', cu.format_number(rmse))]);
            right_element = repmat(' ', 1, 38);
            lines = [left_element '    ' right_element];
            % second row: mae and mape
            left_element = sprintf(['mae:' sprintf('%+34s', cu.format_number(mae))]);
            right_element = sprintf(['mape:' sprintf('%+33s', cu.format_number(mape))]);
            lines = [lines; [left_element '    ' right_element]];
            % third row: Theil's U and bias
            left_element = sprintf(['Theil''s U:' sprintf('%+28s', cu.format_number(theil_u))]);
            right_element = sprintf(['bias:' sprintf('%+33s', cu.format_number(bias))]);
            lines = [lines; [left_element  '    ' right_element]]; 
            % fourth row: log score and crps
            if log_score
                left_element = sprintf(['log score:' sprintf('%+28s', cu.format_number(log_score))]);
                right_element = sprintf(['crps:' sprintf('%+33s', cu.format_number(crps))]);
                lines = [lines; [left_element  '    ' right_element]];
            end
        end
        
        
        function [string] = bool_to_string(bool)

            % [string] = bool_to_string(bool)
            % converts boolean to 'yes' or 'no'
            %
            % parameters  
            % bool: bool
            %     boolean to convert
            %
            % returns
            % string: char
            %     'yes' or 'no', depending on convered boolean  

            if bool
                string = 'yes';
            else
                string = 'no';
            end
        end
        
    end
    
end
        