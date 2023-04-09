classdef Graphics < handle
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------

    
    methods (Access = public) 

        
        function self = Graphics()
        end  
        
        
    end
    

    methods (Access = protected, Hidden = true)
        
   
        function delete_graphics_folder(self)
            % create path to graphics folder
            graphics_folder_path = fullfile(self.project_path, 'graphics');
            % if graphics folder already exists, delete it
            if exist(graphics_folder_path) == 7
                rmdir(graphics_folder_path, 's');
            end
            % create results folder
            mkdir(graphics_folder_path);
        end   
        
        
    end
    
end
