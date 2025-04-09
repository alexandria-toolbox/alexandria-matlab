classdef Tab2RegressionInterface < handle

    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 2 properties (linear regression)
        t2_lr_txt1
        t2_lr_frm1
        t2_lr_bgr1
        t2_lr_rdb1
        t2_lr_rdb2
        t2_lr_rdb3
        t2_lr_rdb4
        t2_lr_rdb5
        t2_lr_rdb6
        t2_lr_txt2
        t2_lr_frm2
        t2_lr_txt3
        t2_lr_edt1
        t2_lr_txt4
        t2_lr_edt2
        t2_lr_txt5
        t2_lr_edt3
        t2_lr_txt6
        t2_lr_frm3
        t2_lr_txt7
        t2_lr_txt8
        t2_lr_edt4
        t2_lr_txt9
        t2_lr_edt5
        t2_lr_txt10
        t2_lr_txt11
        t2_lr_edt6
        t2_lr_txt12
        t2_lr_edt7
        t2_lr_txt13
        t2_lr_txt14
        t2_lr_edt8
        t2_lr_txt15
        t2_lr_edt9
        t2_lr_txt16
        t2_lr_edt10
        t2_lr_txt17
        t2_lr_cbx1
        t2_lr_txt18
        t2_lr_edt11
        t2_lr_txt19
        t2_lr_edt12
        t2_lr_txt20
        t2_lr_txt21
        t2_lr_edt13
        t2_lr_txt22
        t2_lr_edt14
        t2_lr_txt23
        t2_lr_edt15
        t2_lr_txt24
        t2_lr_frm4
        t2_lr_txt25
        t2_lr_cbx2
        t2_lr_txt26
        t2_lr_edt16
        t2_lr_txt27
        t2_lr_edt17
        t2_lr_txt28
        t2_lr_cbx3
        t2_lr_txt29
        t2_lr_edt18
        t2_lr_txt30
        t2_lr_edt19
        t2_lr_txt31
        t2_lr_cbx4
        t2_lr_txt32
        t2_lr_edt20
        t2_lr_txt33
        t2_lr_edt21
        t2_lr_txt34
        t2_lr_frm5
        t2_lr_txt35
        t2_lr_cbx5
        t2_lr_txt36
        t2_lr_cbx6
        t2_lr_txt37
        t2_lr_cbx7
        t2_lr_txt38
        t2_lr_bgr2
        t2_lr_rdb7
        t2_lr_rdb8       
    end
    

    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    
    
    
        function self = Tab2RegressionInterface()
        end
        
        
        function create_tab_2_lr(self)
    
            % regression label
            self.t2_lr_txt1 = uicontrol('style', 'text');
            set(self.t2_lr_txt1, 'unit', 'pixels', 'position', [30 560 300 30]);
            set(self.t2_lr_txt1, 'String', ' Regression type');
            set(self.t2_lr_txt1, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt1, 'FontName', 'Serif');
            set(self.t2_lr_txt1, 'FontSize', 16);
            set(self.t2_lr_txt1, 'FontWeight', 'bold');
            set(self.t2_lr_txt1, 'FontAngle', 'italic');
            set(self.t2_lr_txt1, 'BackgroundColor', self.background_color);  
            set(self.t2_lr_txt1, 'Visible', 'off');

            % frame around regression
            self.t2_lr_frm1 = uicontrol('style','frame');
            set(self.t2_lr_frm1, 'unit', 'pixels', 'position', [20 450 470 110]);
            set(self.t2_lr_frm1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_lr_frm1, 'BackgroundColor', self.background_color);
            set(self.t2_lr_frm1, 'Visible', 'off');

            % regression radiobuttons
            self.t2_lr_bgr1 = uibuttongroup('unit','pixels', 'Position',[25 460 420 90]);
            set(self.t2_lr_bgr1, 'BorderType', 'none');
            set(self.t2_lr_bgr1, 'BackgroundColor', self.background_color); 
            self.t2_lr_rdb1 = uicontrol(self.t2_lr_bgr1,'Style','radiobutton');
            set(self.t2_lr_rdb1, 'Position',[5 70 220 25]);
            set(self.t2_lr_rdb1, 'String',' maximum likelihood');
            set(self.t2_lr_rdb1, 'FontName', 'Serif');
            set(self.t2_lr_rdb1, 'FontSize', 12);
            set(self.t2_lr_rdb1, 'FontWeight', 'bold');
            set(self.t2_lr_rdb1, 'BackgroundColor', self.background_color);
            self.t2_lr_rdb2 = uicontrol(self.t2_lr_bgr1,'Style','radiobutton');
            set(self.t2_lr_rdb2, 'Position',[235 70 220 25]);
            set(self.t2_lr_rdb2, 'String',' simple Bayesian');
            set(self.t2_lr_rdb2, 'FontName', 'Serif');
            set(self.t2_lr_rdb2, 'FontSize', 12);
            set(self.t2_lr_rdb2, 'FontWeight', 'bold');
            set(self.t2_lr_rdb2, 'BackgroundColor', self.background_color);
            self.t2_lr_rdb3 = uicontrol(self.t2_lr_bgr1,'Style','radiobutton');
            set(self.t2_lr_rdb3, 'Position',[5 35 220 25]);
            set(self.t2_lr_rdb3, 'String',' hierarchical');
            set(self.t2_lr_rdb3, 'FontName', 'Serif');
            set(self.t2_lr_rdb3, 'FontSize', 12);
            set(self.t2_lr_rdb3, 'FontWeight', 'bold');
            set(self.t2_lr_rdb3, 'BackgroundColor', self.background_color);
            self.t2_lr_rdb4 = uicontrol(self.t2_lr_bgr1,'Style','radiobutton');
            set(self.t2_lr_rdb4, 'Position',[235 35 220 25]);
            set(self.t2_lr_rdb4, 'String',' independent');
            set(self.t2_lr_rdb4, 'FontName', 'Serif');
            set(self.t2_lr_rdb4, 'FontSize', 12);
            set(self.t2_lr_rdb4, 'FontWeight', 'bold');
            set(self.t2_lr_rdb4, 'BackgroundColor', self.background_color);
            self.t2_lr_rdb5 = uicontrol(self.t2_lr_bgr1,'Style','radiobutton');
            set(self.t2_lr_rdb5, 'Position',[5 0 220 25]);
            set(self.t2_lr_rdb5, 'String',' heteroscedastic');
            set(self.t2_lr_rdb5, 'FontName', 'Serif');
            set(self.t2_lr_rdb5, 'FontSize', 12);
            set(self.t2_lr_rdb5, 'FontWeight', 'bold');
            set(self.t2_lr_rdb5, 'BackgroundColor', self.background_color);
            self.t2_lr_rdb6 = uicontrol(self.t2_lr_bgr1,'Style','radiobutton');
            set(self.t2_lr_rdb6, 'Position',[235 0 220 25]);
            set(self.t2_lr_rdb6, 'String',' autocorrelated');
            set(self.t2_lr_rdb6, 'FontName', 'Serif');
            set(self.t2_lr_rdb6, 'FontSize', 12);
            set(self.t2_lr_rdb6, 'FontWeight', 'bold');
            set(self.t2_lr_rdb6, 'BackgroundColor', self.background_color);
            set(self.t2_lr_bgr1, 'Visible', 'off');
            if self.user_inputs.tab_2_lr.regression_type == 1
                set(self.t2_lr_bgr1, 'SelectedObject', self.t2_lr_rdb1);
            elseif self.user_inputs.tab_2_lr.regression_type == 2
                set(self.t2_lr_bgr1, 'SelectedObject', self.t2_lr_rdb2);
            elseif self.user_inputs.tab_2_lr.regression_type == 3
                set(self.t2_lr_bgr1, 'SelectedObject', self.t2_lr_rdb3);
            elseif self.user_inputs.tab_2_lr.regression_type == 4
                set(self.t2_lr_bgr1, 'SelectedObject', self.t2_lr_rdb4);
            elseif self.user_inputs.tab_2_lr.regression_type == 5
                set(self.t2_lr_bgr1, 'SelectedObject', self.t2_lr_rdb5);
            elseif self.user_inputs.tab_2_lr.regression_type == 6
                set(self.t2_lr_bgr1, 'SelectedObject', self.t2_lr_rdb6);    
            end
            set(self.t2_lr_bgr1, 'SelectionChangeFcn', @self.cb_t2_lr_bgr1);

            % estimation label
            self.t2_lr_txt2 = uicontrol('style', 'text');
            set(self.t2_lr_txt2, 'unit', 'pixels', 'position', [520 560 300 30]);
            set(self.t2_lr_txt2, 'String', ' Estimation');
            set(self.t2_lr_txt2, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt2, 'FontName', 'Serif');
            set(self.t2_lr_txt2, 'FontSize', 16);
            set(self.t2_lr_txt2, 'FontWeight', 'bold');
            set(self.t2_lr_txt2, 'FontAngle', 'italic');
            set(self.t2_lr_txt2, 'BackgroundColor', self.background_color);  
            set(self.t2_lr_txt2, 'Visible', 'off');

            % frame around estimation
            self.t2_lr_frm2 = uicontrol('style','frame');
            set(self.t2_lr_frm2, 'unit', 'pixels', 'position', [510 450 470 110]);
            set(self.t2_lr_frm2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_lr_frm2, 'BackgroundColor', self.background_color);
            set(self.t2_lr_frm2, 'Visible', 'off');

            % iteration label
            self.t2_lr_txt3 = uicontrol('style', 'text');
            set(self.t2_lr_txt3, 'unit', 'pixels', 'position', [520 530 200 20]);
            set(self.t2_lr_txt3, 'String', ' iterations');
            set(self.t2_lr_txt3, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt3, 'FontName', 'Serif');
            set(self.t2_lr_txt3, 'FontSize', 12);
            set(self.t2_lr_txt3, 'FontWeight', 'bold');
            set(self.t2_lr_txt3, 'BackgroundColor', self.background_color); 
            set(self.t2_lr_txt3, 'Visible', 'off');

            % iteration edit
            self.t2_lr_edt1 = uicontrol('style','edit');
            set(self.t2_lr_edt1, 'unit', 'pixels', 'position', [770 527.5 70 25]);
            set(self.t2_lr_edt1, 'HorizontalAlignment', 'center');  
            set(self.t2_lr_edt1, 'Visible', 'off');
            set(self.t2_lr_edt1, 'String', self.user_inputs.tab_2_lr.iterations);
            set(self.t2_lr_edt1, 'CallBack', @self.cb_t2_lr_edt1);

            % burn-in label
            self.t2_lr_txt4 = uicontrol('style', 'text');
            set(self.t2_lr_txt4, 'unit', 'pixels', 'position', [520 495 200 20]);
            set(self.t2_lr_txt4, 'String', ' burn-in');
            set(self.t2_lr_txt4, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt4, 'FontName', 'Serif');
            set(self.t2_lr_txt4, 'FontSize', 12);
            set(self.t2_lr_txt4, 'FontWeight', 'bold');
            set(self.t2_lr_txt4, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt4, 'Visible', 'off');

            % burn-in edit
            self.t2_lr_edt2 = uicontrol('style','edit');
            set(self.t2_lr_edt2, 'unit', 'pixels', 'position', [770 492.5 70 25]);
            set(self.t2_lr_edt2, 'HorizontalAlignment', 'center');  
            set(self.t2_lr_edt2, 'Visible', 'off');
            set(self.t2_lr_edt2, 'String', self.user_inputs.tab_2_lr.burnin);
            set(self.t2_lr_edt2, 'CallBack', @self.cb_t2_lr_edt2);

            % credibility label
            self.t2_lr_txt5 = uicontrol('style', 'text');
            set(self.t2_lr_txt5, 'unit', 'pixels', 'position', [520 460 200 20]);
            set(self.t2_lr_txt5, 'String', ' credibility level');
            set(self.t2_lr_txt5, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt5, 'FontName', 'Serif');
            set(self.t2_lr_txt5, 'FontSize', 12);
            set(self.t2_lr_txt5, 'FontWeight', 'bold');
            set(self.t2_lr_txt5, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt5, 'Visible', 'off');

            % credibility edit
            self.t2_lr_edt3 = uicontrol('style','edit');
            set(self.t2_lr_edt3, 'unit', 'pixels', 'position', [770 457.5 70 25]);
            set(self.t2_lr_edt3, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt3, 'Visible', 'off');
            set(self.t2_lr_edt3, 'String', self.user_inputs.tab_2_lr.model_credibility);
            set(self.t2_lr_edt3, 'CallBack', @self.cb_t2_lr_edt3);

            % hyperparameter label
            self.t2_lr_txt6 = uicontrol('style', 'text');
            set(self.t2_lr_txt6, 'unit', 'pixels', 'position', [30 400 300 30]);
            set(self.t2_lr_txt6, 'String', ' Hyperparameters');
            set(self.t2_lr_txt6, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt6, 'FontName', 'Serif');
            set(self.t2_lr_txt6, 'FontSize', 16);
            set(self.t2_lr_txt6, 'FontWeight', 'bold');
            set(self.t2_lr_txt6, 'FontAngle', 'italic');
            set(self.t2_lr_txt6, 'BackgroundColor', self.background_color);  
            set(self.t2_lr_txt6, 'Visible', 'off');

            % frame around hyperparameters
            self.t2_lr_frm3 = uicontrol('style','frame');
            set(self.t2_lr_frm3, 'unit', 'pixels', 'position', [20 20 470 380]);
            set(self.t2_lr_frm3, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_lr_frm3, 'BackgroundColor', self.background_color);
            set(self.t2_lr_frm3, 'Visible', 'off');

            % Bayesian label
            self.t2_lr_txt7 = uicontrol('style', 'text');
            set(self.t2_lr_txt7, 'unit', 'pixels', 'position', [30 365 300 25]);
            set(self.t2_lr_txt7, 'String', ' All Bayesian');
            set(self.t2_lr_txt7, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt7, 'FontName', 'Serif');
            set(self.t2_lr_txt7, 'FontSize', 14);
            set(self.t2_lr_txt7, 'FontAngle', 'italic');
            set(self.t2_lr_txt7, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt7, 'Visible', 'off');

            % prior mean b label
            self.t2_lr_txt8 = uicontrol('style', 'text');
            set(self.t2_lr_txt8, 'unit', 'pixels', 'position', [30 335 20 20]);
            set(self.t2_lr_txt8, 'String', ' b');
            set(self.t2_lr_txt8, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt8, 'FontName', 'Serif');
            set(self.t2_lr_txt8, 'FontSize', 12);
            set(self.t2_lr_txt8, 'FontWeight', 'bold');
            set(self.t2_lr_txt8, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt8, 'Visible', 'off');

            % prior mean b edit
            self.t2_lr_edt4 = uicontrol('style','edit');
            set(self.t2_lr_edt4, 'unit', 'pixels', 'position', [155 332.5 70 25]);
            set(self.t2_lr_edt4, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt4, 'Visible', 'off');
            set(self.t2_lr_edt4, 'String', self.user_inputs.tab_2_lr.b);
            set(self.t2_lr_edt4, 'CallBack', @self.cb_t2_lr_edt4);

            % prior variance V label
            self.t2_lr_txt9 = uicontrol('style', 'text');
            set(self.t2_lr_txt9, 'unit', 'pixels', 'position', [30 302.5 20 20]);
            set(self.t2_lr_txt9, 'String', ' V');
            set(self.t2_lr_txt9, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt9, 'FontName', 'Serif');
            set(self.t2_lr_txt9, 'FontSize', 12);
            set(self.t2_lr_txt9, 'FontWeight', 'bold');
            set(self.t2_lr_txt9, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt9, 'Visible', 'off');

            % prior variance V edit
            self.t2_lr_edt5 = uicontrol('style','edit');
            set(self.t2_lr_edt5, 'unit', 'pixels', 'position', [155 300 70 25]);
            set(self.t2_lr_edt5, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt5, 'Visible', 'off');
            set(self.t2_lr_edt5, 'String', self.user_inputs.tab_2_lr.V);
            set(self.t2_lr_edt5, 'CallBack', @self.cb_t2_lr_edt5);

            % hierarchical label
            self.t2_lr_txt10 = uicontrol('style', 'text');
            set(self.t2_lr_txt10, 'unit', 'pixels', 'position', [275 370 200 20]);
            set(self.t2_lr_txt10, 'String', ' Hierarchical');
            set(self.t2_lr_txt10, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt10, 'FontName', 'Serif');
            set(self.t2_lr_txt10, 'FontSize', 14);
            set(self.t2_lr_txt10, 'FontAngle', 'italic');
            set(self.t2_lr_txt10, 'BackgroundColor', self.background_color); 
            set(self.t2_lr_txt10, 'Visible', 'off');

            % shape alpha label
            self.t2_lr_txt11 = uicontrol('style', 'text');
            set(self.t2_lr_txt11, 'unit', 'pixels', 'position', [275 335 200 25]);
            set(self.t2_lr_txt11, 'String', ' α');
            set(self.t2_lr_txt11, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt11, 'FontName', 'Serif');
            set(self.t2_lr_txt11, 'FontSize', 16);
            set(self.t2_lr_txt11, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt11, 'Visible', 'off');

            % shape alpha edit
            self.t2_lr_edt6 = uicontrol('style','edit');
            set(self.t2_lr_edt6, 'unit', 'pixels', 'position', [400 332.5 70 25]);
            set(self.t2_lr_edt6, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt6, 'Visible', 'off');
            set(self.t2_lr_edt6, 'String', self.user_inputs.tab_2_lr.alpha);
            set(self.t2_lr_edt6, 'CallBack', @self.cb_t2_lr_edt6);

            % scale delta label
            self.t2_lr_txt12 = uicontrol('style', 'text');
            set(self.t2_lr_txt12, 'unit', 'pixels', 'position', [275 302.5 200 25]);
            set(self.t2_lr_txt12, 'String', ' δ');
            set(self.t2_lr_txt12, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt12, 'FontName', 'Times');
            set(self.t2_lr_txt12, 'FontSize', 15);
            set(self.t2_lr_txt12, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt12, 'Visible', 'off');

            % scale delta edit
            self.t2_lr_edt7 = uicontrol('style','edit');
            set(self.t2_lr_edt7, 'unit', 'pixels', 'position', [400 300 70 25]);
            set(self.t2_lr_edt7, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt7, 'Visible', 'off');
            set(self.t2_lr_edt7, 'String', self.user_inputs.tab_2_lr.delta);
            set(self.t2_lr_edt7, 'CallBack', @self.cb_t2_lr_edt7);

            % heteroscedastic label
            self.t2_lr_txt13 = uicontrol('style', 'text');
            set(self.t2_lr_txt13, 'unit', 'pixels', 'position', [30 262.5 200 20]);
            set(self.t2_lr_txt13, 'String', ' Heteroscedastic');
            set(self.t2_lr_txt13, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt13, 'FontName', 'Serif');
            set(self.t2_lr_txt13, 'FontSize', 14);
            set(self.t2_lr_txt13, 'FontAngle', 'italic');
            set(self.t2_lr_txt13, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt13, 'Visible', 'off');

            % heteroscedastic mean g label
            self.t2_lr_txt14 = uicontrol('style', 'text');
            set(self.t2_lr_txt14, 'unit', 'pixels', 'position', [30 227.5 200 20]);
            set(self.t2_lr_txt14, 'String', ' g');
            set(self.t2_lr_txt14, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt14, 'FontName', 'Serif');
            set(self.t2_lr_txt14, 'FontSize', 12);
            set(self.t2_lr_txt14, 'FontWeight', 'bold');
            set(self.t2_lr_txt14, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt14, 'Visible', 'off');

            % heteroscedastic mean g edit
            self.t2_lr_edt8 = uicontrol('style','edit');
            set(self.t2_lr_edt8, 'unit', 'pixels', 'position', [155 225 70 25]);
            set(self.t2_lr_edt8, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt8, 'Visible', 'off');
            set(self.t2_lr_edt8, 'String', self.user_inputs.tab_2_lr.g);
            set(self.t2_lr_edt8, 'CallBack', @self.cb_t2_lr_edt8);

            % heteroscedastic variance Q label
            self.t2_lr_txt15 = uicontrol('style', 'text');
            set(self.t2_lr_txt15, 'unit', 'pixels', 'position', [30 195 200 20]);
            set(self.t2_lr_txt15, 'String', ' Q');
            set(self.t2_lr_txt15, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt15, 'FontName', 'Serif');
            set(self.t2_lr_txt15, 'FontSize', 12);
            set(self.t2_lr_txt15, 'FontWeight', 'bold');
            set(self.t2_lr_txt15, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt15, 'Visible', 'off');

            % heteroscedastic variance Q edit
            self.t2_lr_edt9 = uicontrol('style','edit');
            set(self.t2_lr_edt9, 'unit', 'pixels', 'position', [155 192.5 70 25]);
            set(self.t2_lr_edt9, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt9, 'Visible', 'off');
            set(self.t2_lr_edt9, 'String', self.user_inputs.tab_2_lr.Q);
            set(self.t2_lr_edt9, 'CallBack', @self.cb_t2_lr_edt9);

            % kernel variance tau label
            self.t2_lr_txt16 = uicontrol('style', 'text');
            set(self.t2_lr_txt16, 'unit', 'pixels', 'position', [30 162.5 200 25]);
            set(self.t2_lr_txt16, 'String', ' τ');
            set(self.t2_lr_txt16, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt16, 'FontName', 'Serif');
            set(self.t2_lr_txt16, 'FontSize', 16);
            set(self.t2_lr_txt16, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt16, 'Visible', 'off');

            % kernel variance tau edit
            self.t2_lr_edt10 = uicontrol('style','edit');
            set(self.t2_lr_edt10, 'unit', 'pixels', 'position', [155 160 70 25]);
            set(self.t2_lr_edt10, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt10, 'Visible', 'off');
            set(self.t2_lr_edt10, 'String', self.user_inputs.tab_2_lr.tau);
            set(self.t2_lr_edt10, 'CallBack', @self.cb_t2_lr_edt10);

            % thinning label
            self.t2_lr_txt17 = uicontrol('style', 'text');
            set(self.t2_lr_txt17, 'unit', 'pixels', 'position', [30 125 200 25]);
            set(self.t2_lr_txt17, 'String', ' thinning');
            set(self.t2_lr_txt17, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt17, 'FontName', 'Serif');
            set(self.t2_lr_txt17, 'FontSize', 12);
            set(self.t2_lr_txt17, 'FontWeight', 'bold');
            set(self.t2_lr_txt17, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt17, 'Visible', 'off');

            % thinning checkbox
            self.t2_lr_cbx1 = uicontrol('style', 'checkbox');
            set(self.t2_lr_cbx1, 'unit', 'pixels', 'position', [210 130 20 20]);
            set(self.t2_lr_cbx1, 'BackgroundColor', self.background_color);
            set(self.t2_lr_cbx1, 'Visible', 'off');
            set(self.t2_lr_cbx1, 'Value', self.user_inputs.tab_2_lr.thinning);
            set(self.t2_lr_cbx1, 'CallBack', @self.cb_t2_lr_cbx1);

            % thinning frequency label
            self.t2_lr_txt18 = uicontrol('style', 'text');
            set(self.t2_lr_txt18, 'unit', 'pixels', 'position', [30 92.5 250 25]);
            set(self.t2_lr_txt18, 'String', ' frequency');
            set(self.t2_lr_txt18, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt18, 'FontName', 'Serif');
            set(self.t2_lr_txt18, 'FontSize', 12);
            set(self.t2_lr_txt18, 'FontWeight', 'bold');
            set(self.t2_lr_txt18, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt18, 'Visible', 'off');

            % thinning frequency edit
            self.t2_lr_edt11 = uicontrol('style','edit');
            set(self.t2_lr_edt11, 'unit', 'pixels', 'position', [155 97.5 70 25]);
            set(self.t2_lr_edt11, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt11, 'Visible', 'off');
            set(self.t2_lr_edt11, 'String', self.user_inputs.tab_2_lr.thinning_frequency);
            set(self.t2_lr_edt11, 'CallBack', @self.cb_t2_lr_edt11);

            % heteroscedasticity regressors Z label
            self.t2_lr_txt19 = uicontrol('style', 'text');
            set(self.t2_lr_txt19, 'unit', 'pixels', 'position', [30 55 300 25]);
            set(self.t2_lr_txt19, 'String', ' Z variables');
            set(self.t2_lr_txt19, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt19, 'FontName', 'Serif');
            set(self.t2_lr_txt19, 'FontSize', 12);
            set(self.t2_lr_txt19, 'FontWeight', 'bold');
            set(self.t2_lr_txt19, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt19, 'Visible', 'off');

            % heteroscedasticity regressors Z edit
            self.t2_lr_edt12 = uicontrol('style','edit');
            set(self.t2_lr_edt12, 'unit', 'pixels', 'position', [35 35 190 25]);
            set(self.t2_lr_edt12, 'HorizontalAlignment', 'left'); 
            set(self.t2_lr_edt12, 'Visible', 'off');
            set(self.t2_lr_edt12, 'String', self.user_inputs.tab_2_lr.Z_variables);
            set(self.t2_lr_edt12, 'CallBack', @self.cb_t2_lr_edt12);

            % autocorrelation label
            self.t2_lr_txt20 = uicontrol('style', 'text');
            set(self.t2_lr_txt20, 'unit', 'pixels', 'position', [275 262.5 200 20]);
            set(self.t2_lr_txt20, 'String', ' Autocorrelated');
            set(self.t2_lr_txt20, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt20, 'FontName', 'Serif');
            set(self.t2_lr_txt20, 'FontSize', 14);
            set(self.t2_lr_txt20, 'FontAngle', 'italic');
            set(self.t2_lr_txt20, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt20, 'Visible', 'off');

            % autocorrelation order q label
            self.t2_lr_txt21 = uicontrol('style', 'text');
            set(self.t2_lr_txt21, 'unit', 'pixels', 'position', [275 227.5 200 20]);
            set(self.t2_lr_txt21, 'String', ' q');
            set(self.t2_lr_txt21, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt21, 'FontName', 'Serif');
            set(self.t2_lr_txt21, 'FontSize', 12);
            set(self.t2_lr_txt21, 'FontWeight', 'bold');
            set(self.t2_lr_txt21, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt21, 'Visible', 'off');

            % autocorrelation order q edit
            self.t2_lr_edt13 = uicontrol('style','edit');
            set(self.t2_lr_edt13, 'unit', 'pixels', 'position', [400 225 70 25]);
            set(self.t2_lr_edt13, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt13, 'Visible', 'off');
            set(self.t2_lr_edt13, 'String', self.user_inputs.tab_2_lr.q);
            set(self.t2_lr_edt13, 'CallBack', @self.cb_t2_lr_edt13);

            % autocorrelation mean p label
            self.t2_lr_txt22 = uicontrol('style', 'text');
            set(self.t2_lr_txt22, 'unit', 'pixels', 'position', [275 195 200 20]);
            set(self.t2_lr_txt22, 'String', ' p');
            set(self.t2_lr_txt22, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt22, 'FontName', 'Serif');
            set(self.t2_lr_txt22, 'FontSize', 12);
            set(self.t2_lr_txt22, 'FontWeight', 'bold');
            set(self.t2_lr_txt22, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt22, 'Visible', 'off');

            % autocorrelation mean p edit
            self.t2_lr_edt14 = uicontrol('style','edit');
            set(self.t2_lr_edt14, 'unit', 'pixels', 'position', [400 192.5 70 25]);
            set(self.t2_lr_edt14, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt14, 'Visible', 'off');
            set(self.t2_lr_edt14, 'String', self.user_inputs.tab_2_lr.p);
            set(self.t2_lr_edt14, 'CallBack', @self.cb_t2_lr_edt14);

            % autocorrelation variance H label
            self.t2_lr_txt23 = uicontrol('style', 'text');
            set(self.t2_lr_txt23, 'unit', 'pixels', 'position', [275 162.5 200 20]);
            set(self.t2_lr_txt23, 'String', ' H');
            set(self.t2_lr_txt23, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt23, 'FontName', 'Serif');
            set(self.t2_lr_txt23, 'FontSize', 12);
            set(self.t2_lr_txt23, 'FontWeight', 'bold');
            set(self.t2_lr_txt23, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt23, 'Visible', 'off');

            % autocorrelation variance H edit
            self.t2_lr_edt15 = uicontrol('style','edit');
            set(self.t2_lr_edt15, 'unit', 'pixels', 'position', [400 160 70 25]);
            set(self.t2_lr_edt15, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt15, 'Visible', 'off');
            set(self.t2_lr_edt15, 'String', self.user_inputs.tab_2_lr.H);
            set(self.t2_lr_edt15, 'CallBack', @self.cb_t2_lr_edt15);

            % exogenous label
            self.t2_lr_txt24 = uicontrol('style', 'text');
            set(self.t2_lr_txt24, 'unit', 'pixels', 'position', [520 400 300 30]);
            set(self.t2_lr_txt24, 'String', ' Exogenous');
            set(self.t2_lr_txt24, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt24, 'FontName', 'Serif');
            set(self.t2_lr_txt24, 'FontSize', 16);
            set(self.t2_lr_txt24, 'FontWeight', 'bold');
            set(self.t2_lr_txt24, 'FontAngle', 'italic');
            set(self.t2_lr_txt24, 'BackgroundColor', self.background_color); 
            set(self.t2_lr_txt24, 'Visible', 'off');

            % frame around exogenous
            self.t2_lr_frm4 = uicontrol('style','frame');
            set(self.t2_lr_frm4, 'unit', 'pixels', 'position', [510 290 470 110]);
            set(self.t2_lr_frm4, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_lr_frm4, 'BackgroundColor', self.background_color);
            set(self.t2_lr_frm4, 'Visible', 'off');

            % constant label
            self.t2_lr_txt25 = uicontrol('style', 'text');
            set(self.t2_lr_txt25, 'unit', 'pixels', 'position', [520 370 200 20]);
            set(self.t2_lr_txt25, 'String', ' constant');
            set(self.t2_lr_txt25, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt25, 'FontName', 'Serif');
            set(self.t2_lr_txt25, 'FontSize', 12);
            set(self.t2_lr_txt25, 'FontWeight', 'bold');
            set(self.t2_lr_txt25, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt25, 'Visible', 'off');

            % constant checkbox
            self.t2_lr_cbx2 = uicontrol('style', 'checkbox');
            set(self.t2_lr_cbx2, 'unit', 'pixels', 'position', [700 370 20 20]);
            set(self.t2_lr_cbx2, 'BackgroundColor', self.background_color);
            set(self.t2_lr_cbx2, 'Visible', 'off');
            set(self.t2_lr_cbx2, 'Value', self.user_inputs.tab_2_lr.constant);
            set(self.t2_lr_cbx2, 'CallBack', @self.cb_t2_lr_cbx2);

            % constant mean label
            self.t2_lr_txt26 = uicontrol('style', 'text');
            set(self.t2_lr_txt26, 'unit', 'pixels', 'position', [740 370 20 20]);
            set(self.t2_lr_txt26, 'String', ' b');
            set(self.t2_lr_txt26, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt26, 'FontName', 'Serif');
            set(self.t2_lr_txt26, 'FontSize', 12);
            set(self.t2_lr_txt26, 'FontWeight', 'bold');
            set(self.t2_lr_txt26, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt26, 'Visible', 'off');

            % constant mean edit
            self.t2_lr_edt16 = uicontrol('style','edit');
            set(self.t2_lr_edt16, 'unit', 'pixels', 'position', [770 367.5 70 25]);
            set(self.t2_lr_edt16, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt16, 'Visible', 'off');
            set(self.t2_lr_edt16, 'String', self.user_inputs.tab_2_lr.b_constant);
            set(self.t2_lr_edt16, 'CallBack', @self.cb_t2_lr_edt16);

            % constant variance label
            self.t2_lr_txt27 = uicontrol('style', 'text');
            set(self.t2_lr_txt27, 'unit', 'pixels', 'position', [870 370 20 20]);
            set(self.t2_lr_txt27, 'String', ' V');
            set(self.t2_lr_txt27, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt27, 'FontName', 'Serif');
            set(self.t2_lr_txt27, 'FontSize', 12);
            set(self.t2_lr_txt27, 'FontWeight', 'bold');
            set(self.t2_lr_txt27, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt27, 'Visible', 'off');

            % constant variance edit
            self.t2_lr_edt17 = uicontrol('style','edit');
            set(self.t2_lr_edt17, 'unit', 'pixels', 'position', [900 367.5 70 25]);
            set(self.t2_lr_edt17, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt17, 'Visible', 'off');
            set(self.t2_lr_edt17, 'String', self.user_inputs.tab_2_lr.V_constant);
            set(self.t2_lr_edt17, 'CallBack', @self.cb_t2_lr_edt17);

            % trend label
            self.t2_lr_txt28 = uicontrol('style', 'text');
            set(self.t2_lr_txt28, 'unit', 'pixels', 'position', [520 335 200 20]);
            set(self.t2_lr_txt28, 'String', ' linear trend');
            set(self.t2_lr_txt28, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt28, 'FontName', 'Serif');
            set(self.t2_lr_txt28, 'FontSize', 12);
            set(self.t2_lr_txt28, 'FontWeight', 'bold');
            set(self.t2_lr_txt28, 'BackgroundColor', self.background_color); 
            set(self.t2_lr_txt28, 'Visible', 'off');

            % trend checkbox
            self.t2_lr_cbx3 = uicontrol('style', 'checkbox');
            set(self.t2_lr_cbx3, 'unit', 'pixels', 'position', [700 335 20 20]);
            set(self.t2_lr_cbx3, 'BackgroundColor', self.background_color);
            set(self.t2_lr_cbx3, 'Visible', 'off');
            set(self.t2_lr_cbx3, 'Value', self.user_inputs.tab_2_lr.trend);
            set(self.t2_lr_cbx3, 'CallBack', @self.cb_t2_lr_cbx3);

            % trend mean label
            self.t2_lr_txt29 = uicontrol('style', 'text');
            set(self.t2_lr_txt29, 'unit', 'pixels', 'position', [740 335 20 20]);
            set(self.t2_lr_txt29, 'String', ' b');
            set(self.t2_lr_txt29, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt29, 'FontName', 'Serif');
            set(self.t2_lr_txt29, 'FontSize', 12);
            set(self.t2_lr_txt29, 'FontWeight', 'bold');
            set(self.t2_lr_txt29, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt29, 'Visible', 'off');

            % trend mean edit
            self.t2_lr_edt18 = uicontrol('style','edit');
            set(self.t2_lr_edt18, 'unit', 'pixels', 'position', [770 332.5 70 25]);
            set(self.t2_lr_edt18, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt18, 'Visible', 'off');
            set(self.t2_lr_edt18, 'String', self.user_inputs.tab_2_lr.b_trend);
            set(self.t2_lr_edt18, 'CallBack', @self.cb_t2_lr_edt18);

            % trend variance label
            self.t2_lr_txt30 = uicontrol('style', 'text');
            set(self.t2_lr_txt30, 'unit', 'pixels', 'position', [870 335 20 20]);
            set(self.t2_lr_txt30, 'String', ' V');
            set(self.t2_lr_txt30, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt30, 'FontName', 'Serif');
            set(self.t2_lr_txt30, 'FontSize', 12);
            set(self.t2_lr_txt30, 'FontWeight', 'bold');
            set(self.t2_lr_txt30, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt30, 'Visible', 'off');

            % trend variance edit
            self.t2_lr_edt19 = uicontrol('style','edit');
            set(self.t2_lr_edt19, 'unit', 'pixels', 'position', [900 332.5 70 25]);
            set(self.t2_lr_edt19, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt19, 'Visible', 'off');
            set(self.t2_lr_edt19, 'String', self.user_inputs.tab_2_lr.V_trend);
            set(self.t2_lr_edt19, 'CallBack', @self.cb_t2_lr_edt19);

            % quadratic trend label
            self.t2_lr_txt31 = uicontrol('style', 'text');
            set(self.t2_lr_txt31, 'unit', 'pixels', 'position', [520 300 200 20]);
            set(self.t2_lr_txt31, 'String', ' quadratic trend');
            set(self.t2_lr_txt31, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt31, 'FontName', 'Serif');
            set(self.t2_lr_txt31, 'FontSize', 12);
            set(self.t2_lr_txt31, 'FontWeight', 'bold');
            set(self.t2_lr_txt31, 'BackgroundColor', self.background_color); 
            set(self.t2_lr_txt31, 'Visible', 'off');

            % quadratic trend checkbox
            self.t2_lr_cbx4 = uicontrol('style', 'checkbox');
            set(self.t2_lr_cbx4, 'unit', 'pixels', 'position', [700 300 20 20]);
            set(self.t2_lr_cbx4, 'BackgroundColor', self.background_color);
            set(self.t2_lr_cbx4, 'Visible', 'off');
            set(self.t2_lr_cbx4, 'Value', self.user_inputs.tab_2_lr.quadratic_trend);
            set(self.t2_lr_cbx4, 'CallBack', @self.cb_t2_lr_cbx4);

            % quadratic trend mean label
            self.t2_lr_txt32 = uicontrol('style', 'text');
            set(self.t2_lr_txt32, 'unit', 'pixels', 'position', [740 300 20 20]);
            set(self.t2_lr_txt32, 'String', ' b');
            set(self.t2_lr_txt32, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt32, 'FontName', 'Serif');
            set(self.t2_lr_txt32, 'FontSize', 12);
            set(self.t2_lr_txt32, 'FontWeight', 'bold');
            set(self.t2_lr_txt32, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt32, 'Visible', 'off');

            % quadratic trend mean edit
            self.t2_lr_edt20 = uicontrol('style','edit');
            set(self.t2_lr_edt20, 'unit', 'pixels', 'position', [770 297.5 70 25]);
            set(self.t2_lr_edt20, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt20, 'Visible', 'off');
            set(self.t2_lr_edt20, 'String', self.user_inputs.tab_2_lr.b_quadratic_trend);
            set(self.t2_lr_edt20, 'CallBack', @self.cb_t2_lr_edt20);

            % quadratic trend variance label
            self.t2_lr_txt33 = uicontrol('style', 'text');
            set(self.t2_lr_txt33, 'unit', 'pixels', 'position', [870 300 20 20]);
            set(self.t2_lr_txt33, 'String', ' V');
            set(self.t2_lr_txt33, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt33, 'FontName', 'Serif');
            set(self.t2_lr_txt33, 'FontSize', 12);
            set(self.t2_lr_txt33, 'FontWeight', 'bold');
            set(self.t2_lr_txt33, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt33, 'Visible', 'off');

            % quadratic trend variance edit
            self.t2_lr_edt21 = uicontrol('style','edit');
            set(self.t2_lr_edt21, 'unit', 'pixels', 'position', [900 297.5 70 25]);
            set(self.t2_lr_edt21, 'HorizontalAlignment', 'center'); 
            set(self.t2_lr_edt21, 'Visible', 'off');
            set(self.t2_lr_edt21, 'String', self.user_inputs.tab_2_lr.V_quadratic_trend);
            set(self.t2_lr_edt21, 'CallBack', @self.cb_t2_lr_edt21);

            % option label
            self.t2_lr_txt34 = uicontrol('style', 'text');
            set(self.t2_lr_txt34, 'unit', 'pixels', 'position', [520 240 300 30]);
            set(self.t2_lr_txt34, 'String', ' Options');
            set(self.t2_lr_txt34, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt34, 'FontName', 'Serif');
            set(self.t2_lr_txt34, 'FontSize', 16);
            set(self.t2_lr_txt34, 'FontWeight', 'bold');
            set(self.t2_lr_txt34, 'FontAngle', 'italic');
            set(self.t2_lr_txt34, 'BackgroundColor', self.background_color); 
            set(self.t2_lr_txt34, 'Visible', 'off');

            % frame around option
            self.t2_lr_frm5 = uicontrol('style','frame');
            set(self.t2_lr_frm5, 'unit', 'pixels', 'position', [510 20 470 220]);
            set(self.t2_lr_frm5, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_lr_frm5, 'BackgroundColor', self.background_color);
            set(self.t2_lr_frm5, 'Visible', 'off');

            % fit label
            self.t2_lr_txt35 = uicontrol('style', 'text');
            set(self.t2_lr_txt35, 'unit', 'pixels', 'position', [520 200 300 20]);
            set(self.t2_lr_txt35, 'String', ' in-sample fit');
            set(self.t2_lr_txt35, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt35, 'FontName', 'Serif');
            set(self.t2_lr_txt35, 'FontSize', 12);
            set(self.t2_lr_txt35, 'FontWeight', 'bold');
            set(self.t2_lr_txt35, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt35, 'Visible', 'off');

            % fit checkbox
            self.t2_lr_cbx5 = uicontrol('style', 'checkbox');
            set(self.t2_lr_cbx5, 'unit', 'pixels', 'position', [930 200 20 20]);
            set(self.t2_lr_cbx5, 'BackgroundColor', self.background_color);
            set(self.t2_lr_cbx5, 'Visible', 'off');
            set(self.t2_lr_cbx5, 'Value', self.user_inputs.tab_2_lr.insample_fit);
            set(self.t2_lr_cbx5, 'CallBack', @self.cb_t2_lr_cbx5);

            % marginal likelihood label
            self.t2_lr_txt36 = uicontrol('style', 'text');
            set(self.t2_lr_txt36, 'unit', 'pixels', 'position', [520 160 300 20]);
            set(self.t2_lr_txt36, 'String', ' marginal likelihood');
            set(self.t2_lr_txt36, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt36, 'FontName', 'Serif');
            set(self.t2_lr_txt36, 'FontSize', 12);
            set(self.t2_lr_txt36, 'FontWeight', 'bold');
            set(self.t2_lr_txt36, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt36, 'Visible', 'off');

            % marginal likelihood checkbox
            self.t2_lr_cbx6 = uicontrol('style', 'checkbox');
            set(self.t2_lr_cbx6, 'unit', 'pixels', 'position', [930 160 20 20]);
            set(self.t2_lr_cbx6, 'BackgroundColor', self.background_color);
            set(self.t2_lr_cbx6, 'Visible', 'off');
            set(self.t2_lr_cbx6, 'Value', self.user_inputs.tab_2_lr.marginal_likelihood);
            set(self.t2_lr_cbx6, 'CallBack', @self.cb_t2_lr_cbx6);

            % optimization label
            self.t2_lr_txt37 = uicontrol('style', 'text');
            set(self.t2_lr_txt37, 'unit', 'pixels', 'position', [520 120 300 20]);
            set(self.t2_lr_txt37, 'String', ' hyperparameter optimization');
            set(self.t2_lr_txt37, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt37, 'FontName', 'Serif');
            set(self.t2_lr_txt37, 'FontSize', 12);
            set(self.t2_lr_txt37, 'FontWeight', 'bold');
            set(self.t2_lr_txt37, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt37, 'Visible', 'off');

            % optimization checkbox
            self.t2_lr_cbx7 = uicontrol('style', 'checkbox');
            set(self.t2_lr_cbx7, 'unit', 'pixels', 'position', [930 120 20 20]);
            set(self.t2_lr_cbx7, 'BackgroundColor', self.background_color);
            set(self.t2_lr_cbx7, 'Visible', 'off');
            set(self.t2_lr_cbx7, 'Value', self.user_inputs.tab_2_lr.hyperparameter_optimization);
            set(self.t2_lr_cbx7, 'CallBack', @self.cb_t2_lr_cbx7);

            % optimization type label
            self.t2_lr_txt38 = uicontrol('style', 'text');
            set(self.t2_lr_txt38, 'unit', 'pixels', 'position', [520 80 300 20]);
            set(self.t2_lr_txt38, 'String', ' optimization type');
            set(self.t2_lr_txt38, 'HorizontalAlignment', 'left');
            set(self.t2_lr_txt38, 'FontName', 'Serif');
            set(self.t2_lr_txt38, 'FontSize', 12);
            set(self.t2_lr_txt38, 'FontWeight', 'bold');
            set(self.t2_lr_txt38, 'BackgroundColor', self.background_color);
            set(self.t2_lr_txt38, 'Visible', 'off');

            % optimization radiobuttons
            self.t2_lr_bgr2 = uibuttongroup('unit','pixels', 'Position',[690 77.5 270 40]);
            set(self.t2_lr_bgr2, 'BorderType', 'none');
            set(self.t2_lr_bgr2, 'BackgroundColor', self.background_color); 
            self.t2_lr_rdb7 = uicontrol(self.t2_lr_bgr2,'Style','radiobutton');
            set(self.t2_lr_rdb7, 'Position',[100 0 220 25]);
            set(self.t2_lr_rdb7, 'String',' simple');
            set(self.t2_lr_rdb7, 'FontName', 'Serif');
            set(self.t2_lr_rdb7, 'FontSize', 12);
            set(self.t2_lr_rdb7, 'FontWeight', 'bold');
            set(self.t2_lr_rdb7, 'BackgroundColor', self.background_color);
            self.t2_lr_rdb8 = uicontrol(self.t2_lr_bgr2,'Style','radiobutton');
            set(self.t2_lr_rdb8, 'Position',[200 0 220 25]);
            set(self.t2_lr_rdb8, 'String',' full');
            set(self.t2_lr_rdb8, 'FontName', 'Serif');
            set(self.t2_lr_rdb8, 'FontSize', 12);
            set(self.t2_lr_rdb8, 'FontWeight', 'bold');
            set(self.t2_lr_rdb8, 'BackgroundColor', self.background_color);
            set(self.t2_lr_bgr2, 'Visible', 'off');
            if self.user_inputs.tab_2_lr.optimization_type == 1
                set(self.t2_lr_bgr2, 'SelectedObject', self.t2_lr_rdb7);
            else
                set(self.t2_lr_bgr2, 'SelectedObject', self.t2_lr_rdb8);
            end
            set(self.t2_lr_bgr2, 'SelectionChangeFcn', @self.cb_t2_lr_bgr2);

            % indicate that tab 2 for linear regression is now created
            self.created_tab_2_lr = true;
        end
        
        
        function hide_tab_2_lr(self)

            % hide all controls            
            set(self.t2_lr_txt1, 'Visible', 'off');
            set(self.t2_lr_frm1, 'Visible', 'off');
            set(self.t2_lr_bgr1, 'Visible', 'off');
            set(self.t2_lr_rdb1, 'Visible', 'off');
            set(self.t2_lr_rdb2, 'Visible', 'off');
            set(self.t2_lr_rdb3, 'Visible', 'off');
            set(self.t2_lr_rdb4, 'Visible', 'off');
            set(self.t2_lr_rdb5, 'Visible', 'off');
            set(self.t2_lr_rdb6, 'Visible', 'off');
            set(self.t2_lr_txt2, 'Visible', 'off');
            set(self.t2_lr_frm2, 'Visible', 'off');
            set(self.t2_lr_txt3, 'Visible', 'off');
            set(self.t2_lr_edt1, 'Visible', 'off');
            set(self.t2_lr_txt4, 'Visible', 'off');
            set(self.t2_lr_edt2, 'Visible', 'off');
            set(self.t2_lr_txt5, 'Visible', 'off');
            set(self.t2_lr_edt3, 'Visible', 'off');
            set(self.t2_lr_txt6, 'Visible', 'off');
            set(self.t2_lr_frm3, 'Visible', 'off');
            set(self.t2_lr_txt7, 'Visible', 'off');
            set(self.t2_lr_txt8, 'Visible', 'off');
            set(self.t2_lr_edt4, 'Visible', 'off');
            set(self.t2_lr_txt9, 'Visible', 'off');
            set(self.t2_lr_edt5, 'Visible', 'off');
            set(self.t2_lr_txt10, 'Visible', 'off');
            set(self.t2_lr_txt11, 'Visible', 'off');
            set(self.t2_lr_edt6, 'Visible', 'off');
            set(self.t2_lr_txt12, 'Visible', 'off');
            set(self.t2_lr_edt7, 'Visible', 'off');
            set(self.t2_lr_txt13, 'Visible', 'off');
            set(self.t2_lr_txt14, 'Visible', 'off');
            set(self.t2_lr_edt8, 'Visible', 'off');
            set(self.t2_lr_txt15, 'Visible', 'off');
            set(self.t2_lr_edt9, 'Visible', 'off');
            set(self.t2_lr_txt16, 'Visible', 'off');
            set(self.t2_lr_edt10, 'Visible', 'off');
            set(self.t2_lr_txt17, 'Visible', 'off');
            set(self.t2_lr_cbx1, 'Visible', 'off');
            set(self.t2_lr_txt18, 'Visible', 'off');
            set(self.t2_lr_edt11, 'Visible', 'off');
            set(self.t2_lr_txt19, 'Visible', 'off');
            set(self.t2_lr_edt12, 'Visible', 'off');
            set(self.t2_lr_txt20, 'Visible', 'off');
            set(self.t2_lr_txt21, 'Visible', 'off');
            set(self.t2_lr_edt13, 'Visible', 'off');
            set(self.t2_lr_txt22, 'Visible', 'off');
            set(self.t2_lr_edt14, 'Visible', 'off');
            set(self.t2_lr_txt23, 'Visible', 'off');
            set(self.t2_lr_edt15, 'Visible', 'off');
            set(self.t2_lr_txt24, 'Visible', 'off');
            set(self.t2_lr_frm4, 'Visible', 'off');
            set(self.t2_lr_txt25, 'Visible', 'off');
            set(self.t2_lr_cbx2, 'Visible', 'off');
            set(self.t2_lr_txt26, 'Visible', 'off');
            set(self.t2_lr_edt16, 'Visible', 'off');
            set(self.t2_lr_txt27, 'Visible', 'off');
            set(self.t2_lr_edt17, 'Visible', 'off');
            set(self.t2_lr_txt28, 'Visible', 'off');
            set(self.t2_lr_cbx3, 'Visible', 'off');
            set(self.t2_lr_txt29, 'Visible', 'off');
            set(self.t2_lr_edt18, 'Visible', 'off');
            set(self.t2_lr_txt30, 'Visible', 'off');
            set(self.t2_lr_edt19, 'Visible', 'off');
            set(self.t2_lr_txt31, 'Visible', 'off');
            set(self.t2_lr_cbx4, 'Visible', 'off');
            set(self.t2_lr_txt32, 'Visible', 'off');
            set(self.t2_lr_edt20, 'Visible', 'off');
            set(self.t2_lr_txt33, 'Visible', 'off');
            set(self.t2_lr_edt21, 'Visible', 'off');
            set(self.t2_lr_txt34, 'Visible', 'off');
            set(self.t2_lr_frm5, 'Visible', 'off');
            set(self.t2_lr_txt35, 'Visible', 'off');
            set(self.t2_lr_cbx5, 'Visible', 'off');
            set(self.t2_lr_txt36, 'Visible', 'off');
            set(self.t2_lr_cbx6, 'Visible', 'off');
            set(self.t2_lr_txt37, 'Visible', 'off');
            set(self.t2_lr_cbx7, 'Visible', 'off');
            set(self.t2_lr_txt38, 'Visible', 'off');
            set(self.t2_lr_bgr2, 'Visible', 'off');
            set(self.t2_lr_rdb7, 'Visible', 'off');
            set(self.t2_lr_rdb8, 'Visible', 'off');

            % update tab color
            set(self.tab_pbt2, 'BackgroundColor', self.backtabs_color);
        end
        
        
        function show_tab_2_lr(self)

            % show all controls
            set(self.t2_lr_txt1, 'Visible', 'on');
            set(self.t2_lr_frm1, 'Visible', 'on');
            set(self.t2_lr_bgr1, 'Visible', 'on');
            set(self.t2_lr_rdb1, 'Visible', 'on');
            set(self.t2_lr_rdb2, 'Visible', 'on');
            set(self.t2_lr_rdb3, 'Visible', 'on');
            set(self.t2_lr_rdb4, 'Visible', 'on');
            set(self.t2_lr_rdb5, 'Visible', 'on');
            set(self.t2_lr_rdb6, 'Visible', 'on');
            set(self.t2_lr_txt2, 'Visible', 'on');
            set(self.t2_lr_frm2, 'Visible', 'on');
            set(self.t2_lr_txt3, 'Visible', 'on');
            set(self.t2_lr_edt1, 'Visible', 'on');
            set(self.t2_lr_txt4, 'Visible', 'on');
            set(self.t2_lr_edt2, 'Visible', 'on');
            set(self.t2_lr_txt5, 'Visible', 'on');
            set(self.t2_lr_edt3, 'Visible', 'on');
            set(self.t2_lr_txt6, 'Visible', 'on');
            set(self.t2_lr_frm3, 'Visible', 'on');
            set(self.t2_lr_txt7, 'Visible', 'on');
            set(self.t2_lr_txt8, 'Visible', 'on');
            set(self.t2_lr_edt4, 'Visible', 'on');
            set(self.t2_lr_txt9, 'Visible', 'on');
            set(self.t2_lr_edt5, 'Visible', 'on');
            set(self.t2_lr_txt10, 'Visible', 'on');
            set(self.t2_lr_txt11, 'Visible', 'on');
            set(self.t2_lr_edt6, 'Visible', 'on');
            set(self.t2_lr_txt12, 'Visible', 'on');
            set(self.t2_lr_edt7, 'Visible', 'on');
            set(self.t2_lr_txt13, 'Visible', 'on');
            set(self.t2_lr_txt14, 'Visible', 'on');
            set(self.t2_lr_edt8, 'Visible', 'on');
            set(self.t2_lr_txt15, 'Visible', 'on');
            set(self.t2_lr_edt9, 'Visible', 'on');
            set(self.t2_lr_txt16, 'Visible', 'on');
            set(self.t2_lr_edt10, 'Visible', 'on');
            set(self.t2_lr_txt17, 'Visible', 'on');
            set(self.t2_lr_cbx1, 'Visible', 'on');
            set(self.t2_lr_txt18, 'Visible', 'on');
            set(self.t2_lr_edt11, 'Visible', 'on');
            set(self.t2_lr_txt19, 'Visible', 'on');
            set(self.t2_lr_edt12, 'Visible', 'on');
            set(self.t2_lr_txt20, 'Visible', 'on');
            set(self.t2_lr_txt21, 'Visible', 'on');
            set(self.t2_lr_edt13, 'Visible', 'on');
            set(self.t2_lr_txt22, 'Visible', 'on');
            set(self.t2_lr_edt14, 'Visible', 'on');
            set(self.t2_lr_txt23, 'Visible', 'on');
            set(self.t2_lr_edt15, 'Visible', 'on');
            set(self.t2_lr_txt24, 'Visible', 'on');
            set(self.t2_lr_frm4, 'Visible', 'on');
            set(self.t2_lr_txt25, 'Visible', 'on');
            set(self.t2_lr_cbx2, 'Visible', 'on');
            set(self.t2_lr_txt26, 'Visible','on');
            set(self.t2_lr_edt16, 'Visible', 'on');
            set(self.t2_lr_txt27, 'Visible', 'on');
            set(self.t2_lr_edt17, 'Visible', 'on');
            set(self.t2_lr_txt28, 'Visible', 'on');
            set(self.t2_lr_cbx3, 'Visible', 'on');
            set(self.t2_lr_txt29, 'Visible', 'on');
            set(self.t2_lr_edt18, 'Visible', 'on');
            set(self.t2_lr_txt30, 'Visible', 'on');
            set(self.t2_lr_edt19, 'Visible', 'on');
            set(self.t2_lr_txt31, 'Visible', 'on');
            set(self.t2_lr_cbx4, 'Visible', 'on');
            set(self.t2_lr_txt32, 'Visible', 'on');
            set(self.t2_lr_edt20, 'Visible', 'on');
            set(self.t2_lr_txt33, 'Visible', 'on');
            set(self.t2_lr_edt21, 'Visible', 'on');
            set(self.t2_lr_txt34, 'Visible', 'on');
            set(self.t2_lr_frm5, 'Visible', 'on');
            set(self.t2_lr_txt35, 'Visible', 'on');
            set(self.t2_lr_cbx5, 'Visible', 'on');
            set(self.t2_lr_txt36, 'Visible', 'on');
            set(self.t2_lr_cbx6, 'Visible', 'on');
            set(self.t2_lr_txt37, 'Visible', 'on');
            set(self.t2_lr_cbx7, 'Visible', 'on');
            set(self.t2_lr_txt38, 'Visible', 'on');
            set(self.t2_lr_bgr2, 'Visible', 'on');
            set(self.t2_lr_rdb7, 'Visible', 'on');
            set(self.t2_lr_rdb8, 'Visible', 'on');
        end
        
        
        function cb_t2_lr_bgr1(self, hObject, callbackdata)
           if get(self.t2_lr_bgr1, 'SelectedObject') == self.t2_lr_rdb1
               self.user_inputs.tab_2_lr.regression_type = 1;
           elseif get(self.t2_lr_bgr1, 'SelectedObject') == self.t2_lr_rdb2
               self.user_inputs.tab_2_lr.regression_type = 2;
           elseif get(self.t2_lr_bgr1, 'SelectedObject') == self.t2_lr_rdb3
               self.user_inputs.tab_2_lr.regression_type = 3;               
           elseif get(self.t2_lr_bgr1, 'SelectedObject') == self.t2_lr_rdb4
               self.user_inputs.tab_2_lr.regression_type = 4;               
           elseif get(self.t2_lr_bgr1, 'SelectedObject') == self.t2_lr_rdb5
               self.user_inputs.tab_2_lr.regression_type = 5;               
           elseif get(self.t2_lr_bgr1, 'SelectedObject') == self.t2_lr_rdb6
               self.user_inputs.tab_2_lr.regression_type = 6;               
           end
        end
        
        
        function cb_t2_lr_edt1(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.iterations = get(self.t2_lr_edt1, 'String');
        end         
        
        
        function cb_t2_lr_edt2(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.burnin = get(self.t2_lr_edt2, 'String');
        end          
        
        
        function cb_t2_lr_edt3(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.model_credibility = get(self.t2_lr_edt3, 'String');
        end            
        
        
        function cb_t2_lr_edt4(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.b = get(self.t2_lr_edt4, 'String');
        end     
        
        
        function cb_t2_lr_edt5(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.V = get(self.t2_lr_edt5, 'String');
        end 

        
        function cb_t2_lr_edt6(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.alpha = get(self.t2_lr_edt6, 'String');
        end         
        
        
        function cb_t2_lr_edt7(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.delta = get(self.t2_lr_edt7, 'String');
        end         
        
        
        function cb_t2_lr_edt8(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.g = get(self.t2_lr_edt8, 'String');
        end          
        
        
        function cb_t2_lr_edt9(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.Q = get(self.t2_lr_edt9, 'String');
        end         
        
        
        function cb_t2_lr_edt10(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.tau = get(self.t2_lr_edt10, 'String');
        end        
        
        
        function cb_t2_lr_cbx1(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.thinning = logical(get(self.t2_lr_cbx1, 'Value'));
        end        
        
        
        function cb_t2_lr_edt11(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.thinning_frequency = get(self.t2_lr_edt11, 'String');
        end          
        
        
        function cb_t2_lr_edt12(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.Z_variables = get(self.t2_lr_edt12, 'String');
        end 
        
        
        function cb_t2_lr_edt13(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.q = get(self.t2_lr_edt13, 'String');
        end         
        
        
        function cb_t2_lr_edt14(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.p = get(self.t2_lr_edt14, 'String');
        end         
        
        
        function cb_t2_lr_edt15(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.H = get(self.t2_lr_edt15, 'String');
        end          
        
        
        function cb_t2_lr_cbx2(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.constant = logical(get(self.t2_lr_cbx2, 'Value'));
        end          
        
        
        function cb_t2_lr_edt16(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.b_constant = get(self.t2_lr_edt16, 'String');
        end        
        
        
        function cb_t2_lr_edt17(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.V_constant = get(self.t2_lr_edt17, 'String');
        end         
        
        
        function cb_t2_lr_cbx3(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.trend = logical(get(self.t2_lr_cbx3, 'Value'));
        end         
        
        
        function cb_t2_lr_edt18(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.b_trend = get(self.t2_lr_edt18, 'String');
        end        
        
        
        function cb_t2_lr_edt19(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.V_trend = get(self.t2_lr_edt19, 'String');
        end         
        
        
        function cb_t2_lr_cbx4(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.quadratic_trend = logical(get(self.t2_lr_cbx4, 'Value'));
        end        
        
 
        function cb_t2_lr_edt20(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.b_quadratic_trend = get(self.t2_lr_edt20, 'String');
        end        
        
        
        function cb_t2_lr_edt21(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.V_quadratic_trend = get(self.t2_lr_edt21, 'String');
        end 
    
        
        function cb_t2_lr_cbx5(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.insample_fit = logical(get(self.t2_lr_cbx5, 'Value'));
        end           
        
        
        function cb_t2_lr_cbx6(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.marginal_likelihood = logical(get(self.t2_lr_cbx6, 'Value'));
        end          
        
        
        function cb_t2_lr_cbx7(self, hObject, callbackdata)
            self.user_inputs.tab_2_lr.hyperparameter_optimization = logical(get(self.t2_lr_cbx7, 'Value'));
        end          
        
        
        function cb_t2_lr_bgr2(self, hObject, callbackdata)
           if get(self.t2_lr_bgr2, 'SelectedObject') == self.t2_lr_rdb7
               self.user_inputs.tab_2_lr.optimization_type = 1;
           elseif get(self.t2_lr_bgr2, 'SelectedObject') == self.t2_lr_rdb8
               self.user_inputs.tab_2_lr.optimization_type = 2;          
           end
        end         
        
    end

end