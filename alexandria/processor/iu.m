classdef iu
    

    % iu stands for input utilities
    % A class containing static methods to handle user inputs such as string, dates, lists, and so on

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function [fixed_char] = fix_char(chr)
    
            % [fixed_char] = fix_char(chr)
            % converts tabs and multiple spaces to single spaces, removes leading and final spaces
            %
            % parameters:
            % char : char
            %     char to fix
            %
            % returns:
            % fixed_char : char
            %     char with spaces made regular

            % convert tabs into spaces
            chr = regexprep(chr, '\t', ' ');
            % convert multiple spaces into single spaces
            chr = regexprep(chr,' +',' ');
            % trim initial and final spaces
            fixed_char = strtrim(chr); 
        end
        
        
        function [string_array] = char_to_array(chr)
            
            % [string_array] = char_to_array(chr)
            % converts char into array of string, with elements split at spaces
            %
            % parameters: 
            % strng: str or array of str
            %     string to split and distribute in array
            %
            % returns:
            % string_array: array of str
            %     array obtained after splitting original char         
            
            % if item is already a string array, ignore and return            
            if isstring(chr)
                string_array = chr;
            % else, item is a char: fix it in case it is incorrectly formated
            else
                chr = iu.fix_char(chr);
                % split it at spaces to convert to string array
                string_array = split(convertCharsToStrings(chr))';
            end
        end
        
        
        function [chr] = array_to_char(string_array)
            
            % [chr] = array_to_char(string_array)
            % converts list of string into char, with elements separated with commas
            %
            % parameters: 
            % string_array: array of str
            %     string array to burst into single char
            %
            % returns:
            % chr: char
            %     char obtained by bursting string array       
            
            % if array is actually a scalar, just convert to string   
            if isa(string_array, 'double') && numel(string_array) == 1
                chr = num2str(string_array);
            elseif isa(string_array, 'double') && numel(string_array) > 1
                chr = regexprep(num2str(string_array'),'\s+',', ');
            % else, item is a string array: convert to single char
            else
                chr = convertStringsToChars(strjoin(string_array,', '));
            end
        end        
        
        
        function [bool] = is_empty(x)
            
            % [bool] = is_empty(x)
            % checks whether input is empty array or empty char/string
            %
            % parameters: 
            % x : array, char or string
            %     elements to be tested for being empty
            %
            % returns:
            % bool: bool
            %     1 if x is empty, 0 otherwise
            
            bool = isempty(x) || (ischar(x) && isempty(x)) || (isstring(x) ...
                   && length(x) == 1 && strlength(x) == 0); 
        end
        
        
        function [bool] = is_digit(x)
            
            % [bool] = is_digit(x)
            % checks whether input is char representing integer value
            % 
            % parameters: 
            % x : char
            %     elements to be tested for being integer value char
            %
            % returns:
            % bool: bool
            %     1 if x is char representing integer, 0 otherwise
            
            if (ischar(x) || isstring(x)) && numel(regexp(x, '\d')) == numel(x)
                bool = true;
            else
                bool = false;
            end
        end
        
        
        function [bool] = is_integer(x)
            
            % [bool] = is_integer(x)
            % checks whether input is double representing integer value
            %
            % parameters:
            % x : double
            %     element to be tested for integer value
            %
            % returns:
            % bool: bool
            %     1 if x is integer, 0 otherwise
            
            if isnumeric(x) && numel(x) == 1 && ~isnan(x) && isfinite(x) && x == floor(x)
                bool = true;
            else
                bool = false;
            end
        end
        
        
        function [bool] = contains_char(string_array, x)
            
            % [bool] = contains_char(string_array, x)
            % checks whether char x appears in string array
            %
            % parameters:
            % string_array : string array
            %     string array potentially including char x
            % x : char
            %     char to be tested for existence in string_array
            %
            % returns:
            % bool : bool
            %     1 if x is in string_array, 0 otherwise
            
            if sum(contains(string_array, x)) >= 1
                bool = true;
            else
                bool = false;
            end
        end
        
        
        function [structure] = concatenate_structures(structure_1, structure_2)
            
            % [structure] = concatenate_structures[structure_1, structure_2]
            % concatenate two structures that have different keys
            %
            % parameters:
            % structure_1 : struct
            %     first structure in concatenation
            % structure_2 : struct
            %     second structure in concatenation  
            %
            % returns:
            % structure : struct
            %     concatenated structure
            
            structure = cell2struct([struct2cell(structure_1);struct2cell(structure_2)], ...
                                    [fieldnames(structure_1);fieldnames(structure_2)]);
        end
            

        function check_file_path(path, file)

            % check_file_path(path, file)
            % checks whether file exists at given path, and is of correct format (csv, xls, xlsx)
            % 
            % parameters: 
            % path : char
            %     path to folder containing data file
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % none

            % check if path is valid
            if ~isfolder(path)
                error(['Path error. Specified path ' path ' does not exist.']);
            end
            % check if file exists
            file_path = fullfile(path, file);
            if ~isfile(file_path)
                error(['File error. File ' file_path ' could not be found.']);
            end
            % check if data file is of correct format (excel or csv)
            file_type = strsplit(file, '.');
            file_type = file_type{end};
            if ~ismember(file_type, ["csv" "xls" "xlsx"])
                error(['Type error for file ' file '. File must be of type csv, xls or xlsx.']);
            end
        end
        
        
        function [data] = load_data(path, file)
            
            % [data] = load_data(path, file)
            % loads data file of type csv, xls or xlsx into Matlab table
            % 
            % parameters: 
            % path : char
            %     path to folder containing data file
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % data : table
            %     table containing loaded data
            
            % make sure path and file strings are properly formatted
            path = iu.fix_char(path);
            file = iu.fix_char(file);
            file_path = fullfile(path, file);
            % load file
            data = readtable(file_path);
            % create index, if relevant
            if isequal(data.Properties.VariableNames{1},'Var1')
                data.Properties.RowNames = string(data{:,'Var1'});
                data = removevars(data,{'Var1'});
            end
        end
        
        
        function check_variables(data, file, variables, tag)

            % check_variables(data, file, variables, tag)
            % checks whether specified variables exist in data table
            %
            % parameters: 
            % data : table
            %     table to examine for data
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % variables : string array
            %     array of endogenous variables to check in table
            % tag : char
            %     tag to apply in error message, if any (e.g. 'Endogenous variables')
            %
            % returns:
            %     none
            
            % obtain the array of variables in table
            table_variables = string(data.Properties.VariableNames);
            % check that variables exist in this array
            missing_variables = setdiff(variables, table_variables);
            if ~iu.is_empty(missing_variables)
                error(['Data error for file ' file '. ' tag ' ' char(strjoin(missing_variables, ', ')) ' cannot be found.']);
            end
        end
        
        
        function check_dates(data, file, start_date, end_date)
            
            % check_dates(data, file, start_date, end_date)
            % check whether start and end dates can be found in index of dataframe
            %
            % parameters:
            % data : table
            %     table to examine for start and end dates
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index
            %
            % returns:
            % none
            
            % first obtain the list of dates (as strings)
            dates = string(data.Properties.RowNames);
            % check for start date
            if ~ismember(start_date, dates)
                error(['Date error for file ' file '. Start date ' convertStringsToChars(start_date) ' cannot be found.']);
            end
            if ~ismember(end_date, dates)
                error(['Date error for file ' file '. End date ' convertStringsToChars(end_date) ' cannot be found.']);
            end
        end

        
        function [sample] = fetch_data(data, file, start_date, end_date, variables, tag)
            
            % [sample] = fetch_data(data, file, start_date, end_date, variables, tag)
            % fetches variables from data at given sample dates
            %
            % parameters:
            % data : table
            %     dataframe containing loaded data
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index    
            % variables : string array
            %     array of variables to search in table
            % tag : char
            %     tag to apply in error message, if any (e.g. 'Endogenous variables')
            %
            % returns:
            % sample : matrix
            %     matrix containing fetched data
            
            % recover sample dates
            dates = string(data.Properties.RowNames);
            start_date_index = find(all(ismember(dates, start_date), 2));
            end_date_index = find(all(ismember(dates, end_date), 2));
            dates = dates(start_date_index:end_date_index);
            % if variable array is not empty, recover sample for given variables and dates
            if ~iu.is_empty(variables)
                sample = data(dates, variables);
                types = convertCharsToStrings(varfun(@class,sample,'OutputFormat','cell'));
                % test for non-numerical values (strings), and if any, raise error
                if any(contains(types,'cell'))
                    error(['Data error for file ' file '. ' tag ' contains text entries, which are unhandled.']);                
                % test for NaNs, and if any, raise error
                elseif any(any(isnan(table2array(sample))))
                    error(['Data error for file ' file '. ' tag ' contains NaN entries, which are unhandled.']);
                % else, data is valid: convert to matrix
                else
                    sample = table2array(sample);
                end
            % if no variable, return empty array
            else
                sample = [];
            end
        end        
        
        
        function [date_format] = infer_date_format(frequency, file, start_date, end_date)
            
            % [date_format] = infer_date_format(frequency, file, start_date, end_date)
            % infer date format for given data file, which can be either periods or timestamps
            %
            % parameters:
            % frequency : int
            %     data frequency, as int between 1 and 6
            % file : char
            %     name of data file (with extension csv, xls or xlsx)            
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index
            %
            % returns:
            % date_format : char
            %     date format, either 'periods' or 'timestamps'
            
            if frequency == 1
                date_format = iu.infer_undated(file, start_date, end_date);
            elseif frequency == 2
                date_format = iu.infer_annual(file, start_date, end_date);
            elseif frequency == 3
                date_format = iu.infer_quarterly(file, start_date, end_date);                
            elseif frequency == 4
                date_format = iu.infer_monthly(file, start_date, end_date);                  
            elseif frequency == 5
                date_format = iu.infer_weekly(file, start_date, end_date);  
            elseif frequency == 6
                date_format = iu.infer_daily(file, start_date, end_date);                 
            end
        end
            
            
        function [date_format] = infer_undated(file, start_date, end_date)
            if ~iu.is_digit(start_date) || ~iu.is_digit(end_date)
                error(['Date error for file ' file '. Unrecognized format for cross-sectional/undated sample start and end dates. Should be integers.']);
            else
                date_format = 'periods';
            end
        end
            
 
        function [date_format] = infer_annual(file, start_date, end_date)
            if iu.is_digit(start_date) && iu.is_digit(end_date)
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for annual sample start and end dates. Should be integers (e.g. 1990) or timestamp (e.g. 1990-12-31).']);
            end
        end

  
        function [date_format] = infer_quarterly(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}Q[1-4]')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}Q[1-4]'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for quarterly sample start and end dates. Should be period (e.g. 1990Q1) or timestamp (e.g. 1990-03-31).']);
            end
        end       
        
 
        function [date_format] = infer_monthly(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}M(0?[1-9]|[1][0-2])$')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}M(0?[1-9]|[1][0-2])$'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for monthly sample start and end dates. Should be period (e.g. 1990M1) or timestamp (e.g. 1990-01-31).']);
            end
        end  

        
        function [date_format] = infer_weekly(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}W(0?[1-9]|[1-4][0-9]|[5][0-3])$')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}W(0?[1-9]|[1-4][0-9]|[5][0-3])$'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for weekly sample start and end dates. Should be period (e.g. 1990W1) or timestamp (e.g. 1990-01-05).']);
            end
        end  
        
        
        function [date_format] = infer_daily(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}D(0?0?[1-9]|0?[1-9]\d|[1-2]\d{2}|3[0-5][0-9]|36[0-6])$')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}D(0?0?[1-9]|0?[1-9]\d|[1-2]\d{2}|3[0-5][0-9]|36[0-6])$'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for daily sample start and end dates. Should be period (e.g. 1990D10) or timestamp (e.g. 1990-01-10).']);
            end
        end
        
        
        function [dates] = generate_dates(data, date_format, frequency, file, start_date, end_date)
            
            % [dates] = generate_dates(data, date_format, frequency, file, start_date, end_date)
            % generates date series, under the form of a date array
            %
            % parameters:
            % data : table
            %     table to examine for start and end dates
            % date_format : char
            %     date format, either 'periods' or 'timestamps'        
            % frequency : int
            %     data frequency, as int between 1 and 6
            % file : char
            %     name of data file (with extension csv, xls or xlsx)            
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index
            %
            % returns:
            % dates : date array
            %     array of datetime entries
            
            % recover sample dates
            dates = string(data.Properties.RowNames);
            start_date_index = find(all(ismember(dates, start_date), 2));
            end_date_index = find(all(ismember(dates, end_date), 2));
            dates = dates(start_date_index:end_date_index);
            % convert to dates, depending on frequency
            if frequency == 1
                dates = iu.generate_undated_dates(dates, date_format, file);
            elseif frequency == 2
                dates = iu.generate_annual_dates(dates, date_format, file);
            elseif frequency == 3
                dates = iu.generate_quarterly_dates(dates, date_format, file);  
            elseif frequency == 4
                dates = iu.generate_monthly_dates(dates, date_format, file);
            elseif frequency == 5
                dates = iu.generate_weekly_dates(dates, date_format, file);                 
            elseif frequency == 6
                dates = iu.generate_daily_dates(dates, date_format, file);                 
            end
        end
        
        
        function [dates] = generate_undated_dates(dates, date_format, file)
            if strcmp(date_format, 'periods')
                try
                    dates = str2num(char(dates));
                catch
                    error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
                end
            else
                error(['Date error for file ' file '. Date format should be periods, not timestamps.']);
            end
        end
                    

        function [dates] = generate_annual_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = datetime(str2num(char(dates))+ 1, 1, 0);
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'end','year');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end
        
        
        function [dates] = generate_quarterly_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    months = str2num(dates(:,6:end)) * 3 + 1;
                    dates = datetime(years, months, 0);
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'end','quarter');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end        
        
        
        function [dates] = generate_monthly_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    months = str2num(dates(:,6:end)) + 1;
                    dates = datetime(years, months, 0);
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'end','month');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end         
        

        function [dates] = generate_weekly_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    days = str2num(dates(:,6:end)) * 7;
                    dates = dateshift(datetime(years, 1, days),'dayofweek','Friday');
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'dayofweek','Friday');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end 

        
        function [dates] = generate_daily_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    days = str2num(dates(:,6:end));
                    dates = datetime(years, 1, days);
                elseif strcmp(date_format, 'timestamps')
                    dates = datetime(dates);
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end
        
        
        function [sample_p] = fetch_forecast_data(data, insample_data, ...
                              variables, file, required, periods, tag)
                          
            % [sample_p] = fetch_forecast_data(data, insample_data, ...
            %              variables, file, required, periods, tag)
            % fetches predictor variables from forecast data file
            %
            % parameters:
            % data : table
            %     table containing loaded data
            % insample_data : matrix
            %     matrix containing in-sample counterparts of forecast variables        
            % variables : string array
            %     array of variables to search in table
            % exogenous_variables : string array
            %     array of variables to extract from table  
            % file : char
            %     name of data file (with extension csv, xls or xlsx) 
            % required : bool
            %     if true, checks that variables are provided in file
            % periods : int
            %     number of forecast periods
            % tag : char
            %     tag indicating the varables to be checked
            %
            % returns:
            % sample_p : matrix
            %     matrix containing forecast data
            
            
            % find any missing variable
            data_variables = string(data.Properties.VariableNames);
            missing_endogenous = setdiff(variables, data_variables);            
            % if some variables are missing
            if ~iu.is_empty(missing_endogenous)            
                % if variables are required and no in-sample data is provided, raise error
                if required && iu.is_empty(insample_data)     
                    error(['Data error for file ' file '. ' tag ' ' char(strjoin(missing_variables, ', ')) ' required to conduct forecast (or forecast evaluation), but cannot be found.'])
                % if variables are required but some in-sample data is provided, replicate final in-sample value
                elseif required && ~iu.is_empty(insample_data)
                    sample_p = repmat(insample_data(end,:), periods, 1);
                % if not required, return empty array
                else
                    sample_p = [];
                end
            % if no variable is missing, recover prediction data 
            else
                % if variables is empty, return empty matrix
                if variables == ""
                    sample_p = [];
                % if there are some variables to fetch, get them
                else
                    sample_p = table2array(data(:,variables));
                    % if too few periods of data are provided, raise error
                    if size(data, 1) < periods
                        error(['Data error for file ' file '. Forecasts must be conducted for ' num2str(periods) ' periods, but ' tag ' is provided for fewer periods.'])
                    % else, reduce data to number of periods
                    else
                        sample_p = sample_p(1:periods,:);
                    end
                    % test for NaNs, and if any, raise error
                    if any(isnan(sample_p(:)))
                        error(['Data error for file ' file '. ' tag ' contains NaN entries, which are unhandled.'])
                    end
                end
            end
        end
                
        
        function [forecast_dates] = generate_forecast_dates(end_date, periods, frequency)
            
            % [forecast_dates] = generate_forecast_dates(end_date, periods, frequency)
            % generates date series for forecasts, under the form of a pandas index, for a given number of periods
            %
            % parameters:
            % end_date : char
            %     sample end date to search in dataframe index            
            % periods : int
            %     number of forecast periods            
            % frequency : int
            %     data frequency, as int between 1 and 6
            %
            % returns:
            % forecast_dates : date array
            %     array of datetime entries
            
            % if frequency is undated, create a range from 1 to periods as forecast dates
            if frequency == 1
                forecast_dates = (1:periods)';
            % if frequency is annual, expand sample by years equal to periods
            elseif frequency == 2
                forecast_dates = (dateshift(end_date + calyears(1:periods),'end','year'))';
            % if frequency is quarterly, expand sample by quarters equal to periods            
            elseif frequency == 3
                forecast_dates = (dateshift(end_date + calquarters(1:periods),'end','quarter'))';
            % if frequency is monthly, expand sample by months equal to periods
            elseif frequency == 4
                forecast_dates = (dateshift(end_date + calmonths(1:periods),'end','month'))';
            % if frequency is weekly, expand sample by weeks equal to periods
            elseif frequency == 5
                forecast_dates = (dateshift(end_date + calweeks(1:periods),'dayofweek','Friday'))';            
            % if frequency is daily, expand sample by days equal to periods
            elseif frequency == 6
                forecast_dates = (end_date + caldays(1:periods))';
            end
        end
        
        
        function check_coefficients_table(data, endogenous_variables, exogenous_variables, ...
                 lags, constant, trend, quadratic_trend, file)
        
            % check_coefficients_table(data, endogenous_variables, exogenous_variables, ...
            %                          lags, constant, trend, quadratic_trend, file)
            % checks whether constrained coefficient table is of valid format
            % 
            % parameters:
            % data : table
            %     table containing constrained coefficient information
            % endogenous_variables : string array
            %     array containing the names of endogenous variables
            % exogenous_variables : string array
            %     array containing the names of exogenous variables
            % lags : int
            %     number of lags for the VAR model
            % constant : bool
            %     set to true if a constant is included in the model
            % trend : bool
            %     set to true if a linear trend is included in the model        
            % quadratic_trend : bool
            %     set to true if a quadratic trend is included in the model
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % none
            
            columns = string(data.Properties.VariableNames);
            if ~isequal(columns, ["variable" "responding_to" "lag" "mean" "variance"])
                error(['Data error for file ' file '. Column names don''t match the required pattern.']);
            end
            types = convertCharsToStrings(varfun(@class,data,'OutputFormat','cell'));
            if ~isequal(types(3), "double") || (isequal(types(3), "double") && ~isequal(data.lag,floor(data.lag)))
                error(['Data error for file ' file '. Some entries in column ' ...
                    '''lag'' are not integers.']);
            end
            if ~isequal(types(4), "double")
                error(['Data error for file ' file '. Some entries in column ' ...
                    '''mean'' are not numeric.']);
            end       
            if ~isequal(types(5), "double")
                error(['Data error for file ' file '. Some entries in column ' ...
                    '''variance'' are not numeric.']);
            end              
            automated_variables = ["constant" "trend" "quadratic_trend"];
            variables = [endogenous_variables exogenous_variables automated_variables];
            all_exogenous = [exogenous_variables automated_variables];
            all_lags = -1:lags;
            for i=1:size(data,1)
                variable_value = char(table2array(data(i,["variable"])));
                responding_value = char(table2array(data(i,["responding_to"])));
                lag_value = table2array(data(i,["lag"]));
                mean_value = table2array(data(i,["mean"]));
                var_value = table2array(data(i,["variance"]));
                if ~iu.contains_char(endogenous_variables, variable_value)
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''variable'', row ' num2str(i) ', does not correspond to an endogenous variable.']);
                end
                if ~iu.contains_char(variables,responding_value)
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''responding_to'', row ' num2str(i) ', does not correspond to any of the model variables.']);
                end
                if isequal(responding_value, 'constant') && ~constant
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''responding_to'', row ' num2str(i) ', is ''constant'', but constant is not activated.']);
                end
                if isequal(responding_value, 'trend') && ~trend
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''responding_to'', row ' num2str(i) ', is ''trend'', but trend is not activated.']);
                end
                if isequal(responding_value, 'quadratic_trend') && ~quadratic_trend
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''responding_to'', row ' num2str(i) ', is ''quadratic_trend'', but quadratic trend is not activated.']);
                end
                if ~isnumeric(lag_value) || isnan(lag_value) || isinf(lag_value) 
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''lag'', row ' num2str(i) ' is not numeric.']);    
                end
                if ~iu.contains_char(all_exogenous,responding_value) && ~ismember(lag_value,all_lags)
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''lag'', row ' num2str(i) ' should be an integer in the range of specified lags.']);
                end                
                if isnan(mean_value) || isinf(mean_value) 
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''mean'', row ' num2str(i) ', is NaN of inf.']);
                end
                if isnan(var_value) || isinf(var_value) 
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''variance'', row ' num2str(i) ', is NaN of inf.']);
                end
                if var_value <= 0
                    error(['Data error for file ' file '. Entry for column ' ...
                    '''variance'', row ' num2str(i) ' should be strictly positive.']);
                end
            end
        end
        
        
        function [constrained_coefficients_table] = get_constrained_coefficients_table(data, ...
                                                    endogenous_variables, exogenous_variables)
        
            % function [constrained_coefficients_table] = get_constrained_coefficients_table(data, ...
            %                                             endogenous_variables, exogenous_variables)
            % recover constrained coefficient table in numeric format
            %
            % parameters:
            % data: table
            %     table containing constrained coefficient information
            % endogenous_variables : string array
            %     array containing the names of endogenous variables
            % exogenous_variables : string array
            %     array containing the names of exogenous variables 
            %
            % returns:
            % constrained_coefficients_table : matrix
            %     ndarray containing numeric values only for constrained coefficients prior

            automated_variables = ["constant" "trend" "quadratic_trend"];
            rows = size(data,1);
            temp = zeros(rows,5);
            for i=1:rows
                variable_value = char(table2array(data(i,["variable"])));
                responding_value = char(table2array(data(i,["responding_to"])));
                lag_value = table2array(data(i,["lag"]));
                mean_value = table2array(data(i,["mean"]));
                var_value = table2array(data(i,["variance"]));   
                temp(i,1) = find(strcmp(endogenous_variables, variable_value));
                if isequal(responding_value, 'constant');
                    temp(i,2) = 0.1;
                elseif isequal(responding_value, 'trend');
                    temp(i,2) = 0.2;
                elseif isequal(responding_value, 'quadratic_trend');
                    temp(i,2) = 0.3;
                elseif iu.contains_char(endogenous_variables,responding_value) 
                    temp(i,2) = find(strcmp(endogenous_variables, responding_value));
                elseif iu.contains_char(exogenous_variables,responding_value) 
                    temp(i,2) = -find(strcmp(exogenous_variables, responding_value));
                end
                if ~iu.contains_char(automated_variables,responding_value)
                    temp(i,3) = lag_value;
                end
                temp(i,4) = mean_value;
                temp(i,5) = var_value;
            end
            constrained_coefficients_table = temp;
        end
        
        
        function check_long_run_table(data, endogenous_variables, file)
            
            % check_long_run_table(data, endogenous_variables, file)
            % checks whether long run prior table is of valid format
            % 
            % parameters:
            % data : table
            %     table containing constrained coefficient information
            % endogenous_variables : string array
            %     array containing the names of endogenous variables
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % none
            
            columns = string(data.Properties.VariableNames);
            if ~isequal(columns, endogenous_variables)
                error(['Data error for file ' file '. Column names don''t match the set of endogenous variables.']);
            end
            rows = string(data.Properties.RowNames)';
            if ~isequal(rows, endogenous_variables)
                error(['Data error for file ' file '. Row names don''t match the set of endogenous variables.']);
            end
            types = convertCharsToStrings(varfun(@class,data,'OutputFormat','cell'));
            for i=1:numel(types)
                if ~isequal(types(i),"double")
                    error(['Data error for file ' file '. Some entries in column ' convertStringsToChars(columns(i)) ' are not numeric.']);
                end
            end
            values = table2array(data);
            dimension = size(data,1);
            for i = 1:dimension
                for j = 1:dimension
                    value = values(i,j);
                    if ~isnumeric(value) || isnan(value) || isinf(value)
                        error(['Data error for file ' file '. Entry in row ' num2str(i) ', column ' num2str(j) ' is not numeric, NaN or inf.']);
                    end
                end
            end
        end
        

        function [text_array] = strseq(text, indices)

            % strseq(text, indices)
            % concatenates indices to text and stores in string array
            %
            % parameters:
            % text : char
            %     text to concatenate to indices
            % indices : vector of size (1,n)
            %     indices to concatenate to text            
            % 
            % returns:
            % text_array : string array
            %     array of concatenated text/indices

            text_array = string();
            for i=1:numel(indices)
                text_array(i) = [text num2str(indices(i))];
            end
        end
        
        function check_condition_table(data, endogenous_variables, periods, file)
            
            % check_condition_table(data, endogenous_variables, periods, file)
            % checks whether condition table is of valid format
            % 
            % parameters:
            % data : table
            %     table containing conditional forcast information
            % endogenous_variables : string array
            %     array containing the names of endogenous variables
            % periods : int
            %     number of forecast periods            
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % none
            
            number_endogenous = numel(endogenous_variables);
            shock_list = iu.strseq('shock',[1:number_endogenous]);
            columns = string(data.Properties.VariableNames);
            expected_columns = ["variable" "period" "mean" "variance" shock_list];
            if ~isequal(columns, expected_columns)
                error(['Data error for file ' file '. Column names don''t match the required pattern.']);
            end
            types = convertCharsToStrings(varfun(@class,data,'OutputFormat','cell'));
            if ~isequal(types(2), "double") || (isequal(types(2), "double") && ~isequal(data.period,floor(data.period)))
                error(['Data error for file ' file '. Some entries in column ' ...
                    '''period'' are not integers.']);
            end
            if ~isequal(types(3), "double")
                error(['Data error for file ' file '. Some entries in column ' ...
                    '''mean'' are not numeric.']);
            end       
            if ~isequal(types(4), "double")
                error(['Data error for file ' file '. Some entries in column ' ...
                    '''variance'' are not numeric.']);
            end
            for i=1:number_endogenous
                if ~isequal(types(4+i), "double")
                    error(['Data error for file ' file '. Entry in row 1, column ' ...
                        '''shock' num2str(i) ''' should be 0 or 1.']);
                end
            end
            rows = size(data,1);
            for i=1:rows
                variable = char(table2array(data(i,["variable"])));
                period = table2array(data(i,["period"]));
                mean = table2array(data(i,["mean"]));
                variance = table2array(data(i,["variance"])); 
                if ~iu.contains_char(endogenous_variables, variable)
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column  ''variable'' does not correspond to one of the model endogenous variables.']);
                end
                if ~isnumeric(period)
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column  ''period'' is not numeric, NaN of inf.']);                    
                end
                if period > periods
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column  ''period'' is larger than the number of forecast periods.']);
                end
                if period <= 0
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column  ''period'' should be a positive integer.']);
                end
                if ~isnumeric(mean)
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column  ''mean'' is not numeric, NaN of inf.']);                    
                end         
                if ~isnumeric(variance)
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column  ''variance'' is not numeric, NaN of inf.']);                    
                end  
                if variance < 0
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column  ''variance'' should be non-negative.']);
                end
                if i == 1
                    for j=1:number_endogenous
                        shock = table2array(data(1,['shock' num2str(j)]));
                        if ~ismember(shock, [0 1])
                            error(['Data error for file ' file '. Entry in row 1, column ' ...
                            '''shock' num2str(i) ''' should be 1 or empty.']);
                        end
                    end
                end
            end
        end
        
        
        function [condition_table, shock_table] = get_condition_table(data, endogenous_variables)
            
            % function [condition_table, shock_table] = get_condition_table(data, endogenous_variables)
            % recover condition table in numeric format
            %
            % parameters:
            % data: table
            %     table containing constrained coefficient information
            % endogenous_variables : string array
            %     array containing the names of endogenous variables
            %
            % returns:
            % condition_table : matrix
            %     ndarray containing numeric values for conditions
            % shock_table : matrix
            %     ndarray containing numeric values for shocks
            
            rows = size(data,1);
            condition_table = zeros(rows,4);
            for i = 1:rows
                variable = char(table2array(data(i,["variable"])));
                condition_table(i,1) = find(strcmp(endogenous_variables, variable));
                condition_table(i,2) = table2array(data(i,["period"]));
                condition_table(i,3) = table2array(data(i,["mean"]));
                if table2array(data(i,["variance"])) == 0
                    condition_table(i,4) = 1e-16;
                else
                    condition_table(i,4) = table2array(data(i,["variance"]));
                end
            end
            shocks = table2array(data(1,5:end))';
            shocks(isnan(shocks)) = 0;
            shock_table = shocks;
        end
        
        
        function [raw_dates] = get_raw_sample_dates(path, file, start_date, end_date)
    
            % [raw_dates] = get_raw_sample_dates(path, file, start_date, end_date)
            % get sample dates, in raw format (as in data file, without any convesion to datetime)
            % 
            % parameters:
            % path : char
            %     path to folder containing data file
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index              
            % 
            % returns:
            % raw_dates : date array
            %     array of datetime entries

            data = iu.load_data(path, file);
            dates = string(data.Properties.RowNames);
            start_date_index = find(all(ismember(dates, start_date), 2));
            end_date_index = find(all(ismember(dates, end_date), 2));
            raw_dates = dates(start_date_index:end_date_index);
        end
        
        
        function check_restriction_table(data, raw_dates, endogenous_variables, proxy_variables, var_type, irf_periods, file)
            
            % check_restriction_table(data, raw_dates, endogenous_variables, irf_periods, file)
            % checks whether restriction table is of valid format
            % 
            % parameters:
            % data : table
            %     table containing constrained coefficient information
            % raw_dates : date array
            %     array of datetime entries               
            % endogenous_variables : string array
            %     array containing the names of endogenous variables
            % proxy_variables : string array
            %     array containing the names of proxy variables            
            % var_type : int
            %     type of VAR model
            % irf_periods : int
            %     number of IRF periods      
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % none        
        
            number_endogenous = numel(endogenous_variables);
            number_proxys = numel(proxy_variables);
            shock_list = iu.strseq('shock',[1:number_endogenous]);
            columns = string(data.Properties.VariableNames);
            expected_columns = ["type" "variable" "period" shock_list];
            if ~isequal(columns, expected_columns)
                error(['Data error for file ' file '. Column names don''t match the required pattern.']);
            end
            rows = size(data,1);
            for i = 1:rows
                restriction_type = char(table2array(data(i,["type"])));
                variable = char(table2array(data(i,["variable"])));
                try
                    period = num2str(table2array(data(i,["period"])));
                catch
                    period = char(table2array(data(i,["period"])));
                end
                if ~ismember(restriction_type, ["sign" "zero" "shock" "historical" "covariance"])
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column ''type'' does not correspond to one of the allowed restriction types.']);                      
                end
                if isequal(restriction_type,'sign') || isequal(restriction_type,'zero')
                    if ~iu.is_digit(period)
                        error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                        ', column ''period'' should be an integer.']); 
                    elseif str2num(period) < 0 || str2num(period) > irf_periods
                        error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                        ', column ''period'' should be an integer between 0 and IRF periods.']); 
                    end
                end
                if (isequal(restriction_type,'shock') || isequal(restriction_type,'historical')) && ~ismember(period,raw_dates)
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column ''period'' should be a sample date.']);                    
                end
                if ~isequal(restriction_type,'shock') && ~isequal(restriction_type,'covariance') && ~ismember(variable,endogenous_variables)
                    error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                    ', column ''variable'' does not correspond to one of the model endogenous variables.']);                     
                end
                for j=1:number_endogenous
                    try
                        shock = num2str(table2array(data(i,3+j)));
                    catch
                        shock = char(table2array(data(i,3+j)));
                    end
                    if ~isequal(shock,'NaN') && ~isempty(shock) && isnan(str2double(shock))
                        error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                        ', column ''shock' num2str(j) ''' is not numeric.']);
                    end
                    if ~ismember(shock,["-1" "0" "1"])
                        error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                        ', column ''shock' num2str(j) ''' is not -1, 0 or 1.']);
                    end
                end
                if ~ismember(numel(nonzeros(table2array(data(i,4:end)))),[1 2])
                    error(['Data error for file ' file '. Ill-defined restrictions in row ' num2str(i) ...
                        ', shock columns must contain either 1 or 2 non-zero coefficients.']);
                end
                if isequal(restriction_type,'covariance') && var_type ~= 7
                    error(['Data error for file ' file '. Covariance restriction found in row ' num2str(i) ', but VAR type is not proxy-SVAR.']);
                elseif isequal(restriction_type,'covariance') && var_type == 7
                    for j=1:(number_endogenous-number_proxys)
                        shock = table2array(data(i,3+j));
                        if shock ~= 0
                            error(['Data error for file ' file '. Entry in row ' num2str(i) ...
                                ', column "shock' num2str(j) '" has covariance restriction while not correlated with proxys.']);
                        end
                    end
                end
            end
        end
        
        
        function [restriction_table] = get_restriction_table(data, raw_dates, endogenous_variables, proxy_variables)
            
            % function [condition_table, shock_table] = get_condition_table(data, endogenous_variables)
            % recover condition table in numeric format
            %
            % parameters:
            % data: table
            %     table containing constrained coefficient information
            % raw_dates : date array
            %     array of datetime entries              
            % endogenous_variables : string array
            %     array containing the names of endogenous variables
            % proxy_variables : string array
            %     array containing the names of proxy variables 
            %
            % returns:
            % restriction_table : matrix
            %     ndarray containing numeric values for restrictions
            
            rows = size(data,1);
            columns = size(data,2);
            number_endogenous = numel(endogenous_variables);
            restriction_types = ["zero" "sign" "shock" "historical" "covariance"];
            restriction_table = zeros(rows,columns);
            for i = 1:rows
                restriction_type = char(table2array(data(i,["type"])));
                period = char(num2str(table2array(data(i,["period"]))));
                variable = char(table2array(data(i,["variable"])));
                restriction_table(i,1) = find(strcmp(restriction_types, restriction_type));
                if isequal(restriction_type,'sign') || isequal(restriction_type,'zero') || isequal(restriction_type,'historical')
                    restriction_table(i,2) = find(strcmp(endogenous_variables, variable));
                elseif isequal(restriction_type,'covariance')
                    restriction_table(i,2) = find(strcmp(proxy_variables, variable));
                end                
                if isequal(restriction_type,'sign') || isequal(restriction_type,'zero')
                    restriction_table(i,3) = str2num(period);
                elseif isequal(restriction_type,'shock') || isequal(restriction_type,'historical')
                    restriction_table(i,3) = find(strcmp(raw_dates, period));
                end
                for j=4:3+number_endogenous
                    restriction_table(i,j) = table2array(data(i,j));
                end
            end
        end


        function [model_name model_class model_type] = identify_model(model)
            
            % identify_model(model)
            % get model and model type
            % 
            % parameters:
            % model : class
            %     class from which model must be extracted
            % 
            % returns:
            % model_name : char
            %     model name
            % model_class : int
            %     general model class (linear regression, VAR, ...)
            % model_type : int
            %     specific model type (maximum likelihood regression, simple Bayesian regression, ...)
            
            class_name = class(model);
            if isequal(class_name, 'MaximumLikelihoodRegression')
                model_name = 'Maximum Likelihood Regression';
                model_class = 1;
                model_type = 1;
            elseif isequal(class_name, 'SimpleBayesianRegression')
                model_name = 'Simple Bayesian Regression';
                model_class = 1;
                model_type = 2;
            elseif isequal(class_name, 'HierarchicalBayesianRegression')
                model_name = 'Hierarchical Bayesian Regression';
                model_class = 1;
                model_type = 3;
            elseif isequal(class_name, 'IndependentBayesianRegression')
                model_name = 'Independent Bayesian Regression';
                model_class = 1;
                model_type = 4;
            elseif isequal(class_name, 'HeteroscedasticBayesianRegression')
                model_name = 'Heteroscedastic Bayesian Regression';
                model_class = 1;
                model_type = 5;
            elseif isequal(class_name, 'AutocorrelatedBayesianRegression')
                model_name = 'Autocorrelated Bayesian Regression';
                model_class = 1;
                model_type = 6;
            elseif isequal(class_name, 'MaximumLikelihoodVar')
                model_name = 'Maximum Likelihood Var';
                model_class = 2;
                model_type = 1;      
            elseif isequal(class_name, 'MinnesotaBayesianVar')
                model_name = 'Minnesota Bayesian Var';
                model_class = 2;
                model_type = 2;
            elseif isequal(class_name, 'NormalWishartBayesianVar')
                model_name = 'Normal-Wishart Bayesian Var';
                model_class = 2;
                model_type = 3;
            elseif isequal(class_name, 'IndependentBayesianVar')
                model_name = 'Independent Bayesian Var';
                model_class = 2;
                model_type = 4;
            elseif isequal(class_name, 'DummyObservationBayesianVar')
                model_name = 'Dummy Observation Bayesian Var';
                model_class = 2;
                model_type = 5;
            elseif isequal(class_name, 'LargeBayesianVar')
                model_name = 'Large Bayesian Var';
                model_class = 2;
                model_type = 6;
            elseif isequal(class_name, 'BayesianProxySvar')
                model_name = 'Bayesian Proxy Svar';
                model_class = 2;
                model_type = 7;               
            end
        end
        
    end
           
end
   