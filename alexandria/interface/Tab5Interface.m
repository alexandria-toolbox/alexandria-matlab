classdef Tab5Interface < handle

    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 5 properties
        t5_txt1
        t5_txt2
        t5_txt3
        t5_txt4
        t5_txt5
        t5_txt6
        t5_txt7
        t5_txt8
        t5_txt9
        t5_txt10
        t5_txt11
        t5_img1
    end
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)
        

        function self = Tab5Interface()
        end
        
        
        function create_tab_5(self)

            % main title
            self.t5_txt1 = uicontrol('style', 'text');
            set(self.t5_txt1, 'unit', 'pixels', 'position', [10 520 300 70]);
            set(self.t5_txt1, 'String', ' Alexandria');
            set(self.t5_txt1, 'HorizontalAlignment', 'left');
            set(self.t5_txt1, 'FontName', 'Serif');
            set(self.t5_txt1, 'FontSize', 34);
            set(self.t5_txt1, 'FontWeight', 'bold');
            set(self.t5_txt1, 'FontAngle', 'italic');
            set(self.t5_txt1, 'BackgroundColor', self.background_color); 
            set(self.t5_txt1, 'Visible', 'off');

            % subtitle
            self.t5_txt2 = uicontrol('style', 'text');
            set(self.t5_txt2, 'unit', 'pixels', 'position', [15 490 300 30]);
            set(self.t5_txt2, 'String', ' By Romain Legrand');
            set(self.t5_txt2, 'HorizontalAlignment', 'left');
            set(self.t5_txt2, 'FontName', 'Serif');
            set(self.t5_txt2, 'FontSize', 16);
            set(self.t5_txt2, 'FontWeight', 'bold');
            set(self.t5_txt2, 'FontAngle', 'italic');
            set(self.t5_txt2, 'BackgroundColor', self.background_color);
            set(self.t5_txt2, 'Visible', 'off');

            % copyright label
            self.t5_txt3 = uicontrol('style', 'text');
            set(self.t5_txt3, 'unit', 'pixels', 'position', [15 450 460 30]);
            set(self.t5_txt3, 'String', ' Copyright © Romain Legrand 2021');
            set(self.t5_txt3, 'HorizontalAlignment', 'left');
            set(self.t5_txt3, 'FontName', 'Serif');
            set(self.t5_txt3, 'FontSize', 16);
            set(self.t5_txt3, 'FontWeight', 'bold');
            set(self.t5_txt3, 'FontAngle', 'italic');
            set(self.t5_txt3, 'BackgroundColor', self.background_color);
            set(self.t5_txt3, 'Visible', 'off');

            % information label
            self.t5_txt4 = uicontrol('style', 'text');
            set(self.t5_txt4, 'unit', 'pixels', 'position', [15 350 460 30]);
            set(self.t5_txt4, 'String', ' for more information:');
            set(self.t5_txt4, 'HorizontalAlignment', 'left');
            set(self.t5_txt4, 'FontName', 'Serif');
            set(self.t5_txt4, 'FontSize', 16);
            set(self.t5_txt4, 'FontWeight', 'bold');
            set(self.t5_txt4, 'FontAngle', 'italic');
            set(self.t5_txt4, 'BackgroundColor', self.background_color);
            set(self.t5_txt4, 'Visible', 'off');

            % mail label
            self.t5_txt5 = uicontrol('style', 'text');
            set(self.t5_txt5, 'unit', 'pixels', 'position', [15 310 460 30]);
            set(self.t5_txt5, 'String', ' alexandria.toolbox@gmail.com');
            set(self.t5_txt5, 'HorizontalAlignment', 'left');
            set(self.t5_txt5, 'FontName', 'Serif');
            set(self.t5_txt5, 'FontSize', 16);
            set(self.t5_txt5, 'FontWeight', 'bold');
            set(self.t5_txt5, 'FontAngle', 'italic');
            set(self.t5_txt5, 'BackgroundColor', self.background_color);
            set(self.t5_txt5, 'Visible', 'off');

            % website label
            self.t5_txt6 = uicontrol('style', 'text');
            set(self.t5_txt6, 'unit', 'pixels', 'position', [15 270 460 30]);
            set(self.t5_txt6, 'String', ' alexandria-toolbox.github.io');
            set(self.t5_txt6, 'HorizontalAlignment', 'left');
            set(self.t5_txt6, 'FontName', 'Serif');
            set(self.t5_txt6, 'FontSize', 16);
            set(self.t5_txt6, 'FontWeight', 'bold');
            set(self.t5_txt6, 'FontAngle', 'italic');
            set(self.t5_txt6, 'BackgroundColor', self.background_color);
            set(self.t5_txt6, 'Visible', 'off');

            % disclaimer line 1
            self.t5_txt7 = uicontrol('style', 'text');
            set(self.t5_txt7, 'unit', 'pixels', 'position', [15 170 460 30]);
            set(self.t5_txt7, 'String', ' Use of this software implies acceptance of the ');
            set(self.t5_txt7, 'HorizontalAlignment', 'left');
            set(self.t5_txt7, 'FontName', 'Serif');
            set(self.t5_txt7, 'FontSize', 11);
            set(self.t5_txt7, 'FontAngle', 'italic');
            set(self.t5_txt7, 'BackgroundColor', self.background_color);
            set(self.t5_txt7, 'Visible', 'off');

            % disclaimer line 2
            self.t5_txt8 = uicontrol('style', 'text');
            set(self.t5_txt8, 'unit', 'pixels', 'position', [15 140 460 30]);
            set(self.t5_txt8, 'String', ' End User Licence Agreement (EULA) for Alexandria.');
            set(self.t5_txt8, 'HorizontalAlignment', 'left');
            set(self.t5_txt8, 'FontName', 'Serif');
            set(self.t5_txt8, 'FontSize', 11);
            set(self.t5_txt8, 'FontAngle', 'italic');
            set(self.t5_txt8, 'BackgroundColor', self.background_color);
            set(self.t5_txt8, 'Visible', 'off');

            % disclaimer line 3
            self.t5_txt9 = uicontrol('style', 'text');
            set(self.t5_txt9, 'unit', 'pixels', 'position', [15 110 460 30]);
            set(self.t5_txt9, 'String', ' Please read and accept the End-User License ');
            set(self.t5_txt9, 'HorizontalAlignment', 'left');
            set(self.t5_txt9, 'FontName', 'Serif');
            set(self.t5_txt9, 'FontSize', 11);
            set(self.t5_txt9, 'FontAngle', 'italic');
            set(self.t5_txt9, 'BackgroundColor', self.background_color);
            set(self.t5_txt9, 'Visible', 'off');

            % disclaimer line 4
            self.t5_txt10 = uicontrol('style', 'text');
            set(self.t5_txt10, 'unit', 'pixels', 'position', [15 80 460 30]);
            set(self.t5_txt10, 'String', ' Agreement carefully before downloading, installing');
            set(self.t5_txt10, 'HorizontalAlignment', 'left');
            set(self.t5_txt10, 'FontName', 'Serif');
            set(self.t5_txt10, 'FontSize', 11);
            set(self.t5_txt10, 'FontAngle', 'italic');
            set(self.t5_txt10, 'BackgroundColor', self.background_color);
            set(self.t5_txt10, 'Visible', 'off');

            % disclaimer line 5
            self.t5_txt11 = uicontrol('style', 'text');
            set(self.t5_txt11, 'unit', 'pixels', 'position', [15 50 460 30]);
            set(self.t5_txt11, 'String', ' or using Alexandria.');
            set(self.t5_txt11, 'HorizontalAlignment', 'left');
            set(self.t5_txt11, 'FontName', 'Serif');
            set(self.t5_txt11, 'FontSize', 11);
            set(self.t5_txt11, 'FontAngle', 'italic');
            set(self.t5_txt11, 'BackgroundColor', self.background_color);
            set(self.t5_txt11, 'Visible', 'off');

            % side image
            self.t5_img1 = axes('units', 'pixels', 'position', [480 125 480 360]);
            [img, map, alphachannel] = imread(fullfile(self.interface_path, 'credits.png'));
            image(img);
            set(self.t5_img1, 'xtick', [], 'ytick', []);
            set(self.t5_img1, 'Visible', 'off');
            set(get(self.t5_img1,'children'),'visible','off');
        end
        
        
        function hide_tab_5(self)

            % hide all controls            
            set(self.t5_txt1, 'Visible', 'off');
            set(self.t5_txt2, 'Visible', 'off');
            set(self.t5_txt3, 'Visible', 'off');
            set(self.t5_txt4, 'Visible', 'off');
            set(self.t5_txt5, 'Visible', 'off');
            set(self.t5_txt6, 'Visible', 'off');
            set(self.t5_txt7, 'Visible', 'off');
            set(self.t5_txt8, 'Visible', 'off');
            set(self.t5_txt9, 'Visible', 'off');
            set(self.t5_txt10, 'Visible', 'off');
            set(self.t5_txt11, 'Visible', 'off');
            set(self.t5_img1, 'Visible', 'off');
            set(get(self.t5_img1,'children'),'visible','off');

            % update tab color
            set(self.tab_pbt5, 'BackgroundColor', self.backtabs_color);
        end
        
        
        function show_tab_5(self)
        
            % show all controls
            set(self.t5_txt1, 'Visible', 'on');
            set(self.t5_txt2, 'Visible', 'on');
            set(self.t5_txt3, 'Visible', 'on');
            set(self.t5_txt4, 'Visible', 'on');
            set(self.t5_txt5, 'Visible', 'on');
            set(self.t5_txt6, 'Visible', 'on');
            set(self.t5_txt7, 'Visible', 'on');
            set(self.t5_txt8, 'Visible', 'on');
            set(self.t5_txt9, 'Visible', 'on');
            set(self.t5_txt10, 'Visible', 'on');
            set(self.t5_txt11, 'Visible', 'on');
            set(self.t5_img1, 'Visible', 'on');
            set(get(self.t5_img1, 'children'), 'visible', 'on');
        end
        
        
    end
    
    
end