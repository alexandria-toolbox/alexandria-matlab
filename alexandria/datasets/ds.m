classdef ds
    

    % ds stands for DataSets
    % A class containing static methods to load datasets

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)

        
        function [data] = load_taylor_table()

            % [data] = load_taylor_table()
            % load the Taylor dataset as a Matlab table
            %
            % parameters:
            % none
            %
            % returns:
            % data : Matlab table
            %     table containing the Taylor dataset

            self = ds();
            file_path = self.get_file_path('taylor');
            opts = detectImportOptions(file_path);
            opts.VariableTypes{1} = 'char';
            opts.VariableTypes([2:end]) = {'double'};
            data = readtable(file_path, opts, 'ReadRowNames',true);
            data = removevars(data,{'Var1'});
        end
        
        
        function [data] = load_taylor()

            % load_taylor()
            % load the raw Taylor dataset as a matrix
            %
            % parameters:
            % none
            %
            % returns:
            % data : matrix
            %     matrix containing the raw data for the Taylor dataset

            datatable = ds.load_taylor_table();
            data = table2array(datatable);
        end
        
    end
    
    
    methods (Access = protected, Hidden = true)    

        
        function [file_path] = get_file_path(self, file_name)
            dataset_folder_path = fileparts(which(mfilename));
            file_path = fullfile(dataset_folder_path, [file_name '.csv']);
        end
        
    end
    
end