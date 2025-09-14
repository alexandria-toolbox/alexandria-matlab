classdef DataSets < handle
    

    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
    end    


    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    


        function self = DataSets()          
        end


        function [data] = load_taylor_table(self)

            % [data] = load_taylor_table()
            % load the Taylor dataset as a Matlab table
            %
            % parameters:
            % none
            %
            % returns:
            % data : Matlab table
            %     table containing the Taylor dataset

            file_path = self.get_file_path('taylor');
            opts = detectImportOptions(file_path);
            opts.VariableTypes{1} = 'char';
            opts.VariableTypes([2:end]) = {'double'};
            data = readtable(file_path, opts, 'ReadRowNames',true);
            data = removevars(data,{'Var1'});
        end
        
        
        function [data] = load_taylor(self)

            % load_taylor()
            % load the raw Taylor dataset as a matrix
            %
            % parameters:
            % none
            %
            % returns:
            % data : matrix
            %     matrix containing the raw data for the Taylor dataset

            datatable = self.load_taylor_table();
            data = table2array(datatable);
        end


        function [data] = load_islm_table(self)

            % [data] = load_islm_table()
            % load the Euro Area IS-LM as a Matlab table
            %
            % parameters:
            % none
            %
            % returns:
            % data : Matlab table
            %     table containing the IS-LM dataset

            file_path = self.get_file_path('islm');
            opts = detectImportOptions(file_path);
            opts.VariableTypes{1} = 'char';
            opts.VariableTypes([2:end]) = {'double'};
            data = readtable(file_path, opts, 'ReadRowNames',true);
            data = removevars(data,{'Var1'});
        end
        
        
        function [data] = load_islm(self)

            % load_islm()
            % load the raw Euro Area IS-LM dataset as a matrix
            %
            % parameters:
            % none
            %
            % returns:
            % data : matrix
            %     matrix containing the raw data for the IS-LM dataset

            datatable = self.load_islm_table();
            data = table2array(datatable);
        end


        function [data] = load_fdi_table(self)

            % [data] = load_fdi_table()
            % load the India FDI dataset as a Matlab table
            %
            % parameters:
            % none
            %
            % returns:
            % data : Matlab table
            %     table containing the FDI dataset

            file_path = self.get_file_path('fdi');
            opts = detectImportOptions(file_path);
            opts.VariableTypes{1} = 'char';
            opts.VariableTypes([2:end]) = {'double'};
            data = readtable(file_path, opts, 'ReadRowNames',true);
            data = removevars(data,{'Var1'});
        end
        
        
        function [data] = load_fdi(self)

            % load_fdi()
            % load the raw India FDI dataset as a matrix
            %
            % parameters:
            % none
            %
            % returns:
            % data : matrix
            %     matrix containing the raw data for the FDI dataset

            datatable = self.load_fdi_table();
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