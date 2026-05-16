classdef Tab2NowcastingInterface < handle

    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 2 properties (nowcasting)
        t2_now_txt1
        t2_now_txt2        
        t2_now_txt3
        t2_now_txt4        
        t2_now_txt5        
        t2_now_txt6
        t2_now_txt7
        t2_now_txt8
        t2_now_txt9
        t2_now_txt10
        t2_now_txt11
        t2_now_txt12
        t2_now_txt13
        t2_now_txt14
        t2_now_txt15
        t2_now_txt16
        t2_now_txt17
        t2_now_txt18
        t2_now_txt19
        t2_now_txt20
        t2_now_txt21
        t2_now_txt22
        t2_now_txt23
        t2_now_txt24
        t2_now_txt25
        t2_now_txt26
        t2_now_txt27
        t2_now_txt28
        t2_now_txt29
        t2_now_txt30
        t2_now_txt31
        t2_now_txt32
        t2_now_txt33
        t2_now_txt34
        t2_now_txt35
        t2_now_txt36
        t2_now_txt37
        t2_now_txt38
        t2_now_txt39
        t2_now_frm1
        t2_now_frm2  
        t2_now_frm3  
        t2_now_frm4         
        t2_now_bgr1
        t2_now_rdb1
        t2_now_rdb2
        t2_now_rdb3
        t2_now_edt1
        t2_now_edt2
        t2_now_edt3
        t2_now_edt4
        t2_now_edt5
        t2_now_edt6
        t2_now_edt7
        t2_now_edt8
        t2_now_edt9
        t2_now_edt10
        t2_now_edt11
        t2_now_edt12
        t2_now_edt13
        t2_now_edt14
        t2_now_edt15
        t2_now_edt16
        t2_now_edt17
        t2_now_edt18
        t2_now_edt19
        t2_now_edt20
        t2_now_edt21
        t2_now_edt22
        t2_now_edt23
        t2_now_edt24
        t2_now_edt25
        t2_now_edt26
        t2_now_edt27
        t2_now_edt28
        t2_now_mnu1
        t2_now_cbx1
        t2_now_cbx2
        t2_now_cbx3
        t2_now_cbx4
    end
    

    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    
    
    
        function self = Tab2NowcastingInterface()
        end

        
        function create_tab_2_now(self)

            % model label
            self.t2_now_txt1 = uicontrol('style', 'text');
            set(self.t2_now_txt1, 'unit', 'pixels', 'position', [30 560 300 30]);
            set(self.t2_now_txt1, 'String', ' Model');
            set(self.t2_now_txt1, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt1, 'FontName', 'Serif');
            set(self.t2_now_txt1, 'FontSize', 15);
            set(self.t2_now_txt1, 'FontWeight', 'bold');
            set(self.t2_now_txt1, 'FontAngle', 'italic');
            set(self.t2_now_txt1, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt1, 'Visible', 'off');

            % frame around model
            self.t2_now_frm1 = uicontrol('style','frame');
            set(self.t2_now_frm1, 'unit', 'pixels', 'position', [20 415 470 145]);
            set(self.t2_now_frm1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_now_frm1, 'BackgroundColor', self.background_color);
            set(self.t2_now_frm1, 'Visible', 'off');

            % Selection label
            self.t2_now_txt2 = uicontrol('style', 'text');
            set(self.t2_now_txt2, 'unit', 'pixels', 'position', [30 525 200 30]);
            set(self.t2_now_txt2, 'String', ' Selection');
            set(self.t2_now_txt2, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt2, 'FontName', 'Serif');
            set(self.t2_now_txt2, 'FontSize', 14);
            set(self.t2_now_txt2, 'FontAngle', 'italic');
            set(self.t2_now_txt2, 'BackgroundColor', self.background_color);
            set(self.t2_now_txt2, 'Visible', 'off');  

            % model radiobuttons
            % self.t2_now_bgr1 = uibuttongroup('unit','pixels', 'Position',[25 455 440 75]);
            self.t2_now_bgr1 = uibuttongroup('unit','pixels', 'Position',[25 450 440 80]);
            set(self.t2_now_bgr1, 'BorderType', 'none');
            set(self.t2_now_bgr1, 'BackgroundColor', self.background_color); 
            self.t2_now_rdb1 = uicontrol(self.t2_now_bgr1,'Style','radiobutton');
            set(self.t2_now_rdb1, 'Position',[5 53 220 25]);
            set(self.t2_now_rdb1, 'String',' mixed frequency BVAR');
            set(self.t2_now_rdb1, 'FontName', 'Serif');
            set(self.t2_now_rdb1, 'FontSize', 11);
            set(self.t2_now_rdb1, 'FontWeight', 'bold');
            set(self.t2_now_rdb1, 'BackgroundColor', self.background_color);
            self.t2_now_rdb2 = uicontrol(self.t2_now_bgr1,'Style','radiobutton');
            set(self.t2_now_rdb2, 'Position',[5 27 220 25]);
            set(self.t2_now_rdb2, 'String',' dynamic factor model');
            set(self.t2_now_rdb2, 'FontName', 'Serif');
            set(self.t2_now_rdb2, 'FontSize', 11);
            set(self.t2_now_rdb2, 'FontWeight', 'bold');
            set(self.t2_now_rdb2, 'BackgroundColor', self.background_color);
            self.t2_now_rdb3 = uicontrol(self.t2_now_bgr1,'Style','radiobutton');
            set(self.t2_now_rdb3, 'Position',[5 0 220 25]);
            set(self.t2_now_rdb3, 'String',' Midas regression');
            set(self.t2_now_rdb3, 'FontName', 'Serif');
            set(self.t2_now_rdb3, 'FontSize', 11);
            set(self.t2_now_rdb3, 'FontWeight', 'bold');
            set(self.t2_now_rdb3, 'BackgroundColor', self.background_color);
            set(self.t2_now_bgr1, 'Visible', 'off');
            if self.user_inputs.tab_2_now.model == 1
                set(self.t2_now_bgr1, 'SelectedObject', self.t2_now_rdb1);
            elseif self.user_inputs.tab_2_now.model == 2
                set(self.t2_now_bgr1, 'SelectedObject', self.t2_now_rdb2);
            elseif self.user_inputs.tab_2_now.model == 3
                set(self.t2_now_bgr1, 'SelectedObject', self.t2_now_rdb3);
            end
            set(self.t2_now_bgr1, 'SelectionChangeFcn', @self.cb_t2_now_bgr1);

            % Gibbs sampling label
            self.t2_now_txt3 = uicontrol('style', 'text');
            set(self.t2_now_txt3, 'unit', 'pixels', 'position', [260 525 200 30]);
            set(self.t2_now_txt3, 'String', ' Gibbs sampling');
            set(self.t2_now_txt3, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt3, 'FontName', 'Serif');
            set(self.t2_now_txt3, 'FontSize', 14);
            set(self.t2_now_txt3, 'FontAngle', 'italic');
            set(self.t2_now_txt3, 'BackgroundColor', self.background_color);
            set(self.t2_now_txt3, 'Visible', 'off');

            % iteration label
            self.t2_now_txt4 = uicontrol('style', 'text');
            set(self.t2_now_txt4, 'unit', 'pixels', 'position', [260 498 200 25]);
            set(self.t2_now_txt4, 'String', ' iterations');
            set(self.t2_now_txt4, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt4, 'FontName', 'Serif');
            set(self.t2_now_txt4, 'FontSize', 11);
            set(self.t2_now_txt4, 'FontWeight', 'bold');
            set(self.t2_now_txt4, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt4, 'Visible', 'off');  

            % iteration edit
            self.t2_now_edt1 = uicontrol('style','edit');
            set(self.t2_now_edt1, 'unit', 'pixels', 'position', [410 503 70 23]);
            set(self.t2_now_edt1, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt1, 'Visible', 'off');
            set(self.t2_now_edt1, 'String', self.user_inputs.tab_2_now.iterations);
            set(self.t2_now_edt1, 'CallBack', @self.cb_t2_now_edt1);

            % burn-in label
            self.t2_now_txt5 = uicontrol('style', 'text');
            set(self.t2_now_txt5, 'unit', 'pixels', 'position', [260 471 200 25]);
            set(self.t2_now_txt5, 'String', ' burn-in');
            set(self.t2_now_txt5, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt5, 'FontName', 'Serif');
            set(self.t2_now_txt5, 'FontSize', 11);
            set(self.t2_now_txt5, 'FontWeight', 'bold');
            set(self.t2_now_txt5, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt5, 'Visible', 'off');   

            % burn-in edit
            self.t2_now_edt2 = uicontrol('style','edit');
            set(self.t2_now_edt2, 'unit', 'pixels', 'position', [410 476 70 23]);
            set(self.t2_now_edt2, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt2, 'Visible', 'off');
            set(self.t2_now_edt2, 'String', self.user_inputs.tab_2_now.burnin);
            set(self.t2_now_edt2, 'CallBack', @self.cb_t2_now_edt2);

            % credibility label
            self.t2_now_txt6 = uicontrol('style', 'text');
            set(self.t2_now_txt6, 'unit', 'pixels', 'position', [260 451 200 20]);
            set(self.t2_now_txt6, 'String', ' credibility level');
            set(self.t2_now_txt6, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt6, 'FontName', 'Serif');
            set(self.t2_now_txt6, 'FontSize', 11);
            set(self.t2_now_txt6, 'FontWeight', 'bold');
            set(self.t2_now_txt6, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt6, 'Visible', 'off');         

            % credibility edit
            self.t2_now_edt3 = uicontrol('style','edit');
            set(self.t2_now_edt3, 'unit', 'pixels', 'position', [410 450 70 23]);
            set(self.t2_now_edt3, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt3, 'Visible', 'off');
            set(self.t2_now_edt3, 'String', self.user_inputs.tab_2_now.model_credibility);
            set(self.t2_now_edt3, 'CallBack', @self.cb_t2_now_edt3);

            % midas label
            self.t2_now_txt7 = uicontrol('style', 'text');
            set(self.t2_now_txt7, 'unit', 'pixels', 'position', [520 560 400 30]);
            set(self.t2_now_txt7, 'String', ' Bayesian MIDAS regression');
            set(self.t2_now_txt7, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt7, 'FontName', 'Serif');
            set(self.t2_now_txt7, 'FontSize', 15);
            set(self.t2_now_txt7, 'FontWeight', 'bold');
            set(self.t2_now_txt7, 'FontAngle', 'italic');
            set(self.t2_now_txt7, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt7, 'Visible', 'off');
            
            % frame around midas regression
            self.t2_now_frm2 = uicontrol('style','frame');
            set(self.t2_now_frm2, 'unit', 'pixels', 'position', [510 415 470 145]);
            set(self.t2_now_frm2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_now_frm2, 'BackgroundColor', self.background_color);
            set(self.t2_now_frm2, 'Visible', 'off');

            % prior label
            self.t2_now_txt11 = uicontrol('style', 'text');
            set(self.t2_now_txt11, 'unit', 'pixels', 'position', [520 530 200 20]);
            set(self.t2_now_txt11, 'String', ' model');
            set(self.t2_now_txt11, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt11, 'FontName', 'Serif');
            set(self.t2_now_txt11, 'FontSize', 11);
            set(self.t2_now_txt11, 'FontWeight', 'bold');
            set(self.t2_now_txt11, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt11, 'Visible', 'off'); 

            % model selection menu
            self.t2_now_mnu1 = uicontrol('style', 'popupmenu');
            set(self.t2_now_mnu1, 'position',[680 530 290 23]);            
            set(self.t2_now_mnu1, 'String', {'unrestricted - minnesota', 'unrestricted - horseshoe', ...
                                              'unrestricted - lasso', 'almon - minnesota', ...
                                              'almon - horseshoe', 'almon - lasso', ...
                                              'fourier - minnesota', 'fourier - horseshoe', ...
                                              'fourier - lasso'});  
            set(self.t2_now_mnu1, 'Visible', 'off');
            % set(self.t2_now_mnu1, 'Value', self.user_inputs.tab_2_now.midas_prior_type);
            set(self.t2_now_mnu1, 'Value', self.user_inputs.tab_2_now.midas_model);
            set(self.t2_now_mnu1, 'CallBack', @self.cb_t2_now_mnu1);

            % endogenous lags label
            self.t2_now_txt8 = uicontrol('style', 'text');
            set(self.t2_now_txt8, 'unit', 'pixels', 'position', [520 498 200 25]);
            set(self.t2_now_txt8, 'String', ' endogenous lags');
            set(self.t2_now_txt8, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt8, 'FontName', 'Serif');
            set(self.t2_now_txt8, 'FontSize', 11);
            set(self.t2_now_txt8, 'FontWeight', 'bold');
            set(self.t2_now_txt8, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt8, 'Visible', 'off');  

            % endogenous lags edit
            self.t2_now_edt4 = uicontrol('style','edit');
            set(self.t2_now_edt4, 'unit', 'pixels', 'position', [680 503 70 23]);
            set(self.t2_now_edt4, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt4, 'Visible', 'off');
            set(self.t2_now_edt4, 'String', self.user_inputs.tab_2_now.midas_endogenous_lags);
            set(self.t2_now_edt4, 'CallBack', @self.cb_t2_now_edt4);

            % exogenous lags label
            self.t2_now_txt9 = uicontrol('style', 'text');
            set(self.t2_now_txt9, 'unit', 'pixels', 'position', [750 498 200 25]);
            set(self.t2_now_txt9, 'String', ' exogenous lags');
            set(self.t2_now_txt9, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt9, 'FontName', 'Serif');
            set(self.t2_now_txt9, 'FontSize', 11);
            set(self.t2_now_txt9, 'FontWeight', 'bold');
            set(self.t2_now_txt9, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt9, 'Visible', 'off');  

            % exogenous lags edit
            self.t2_now_edt5 = uicontrol('style','edit');
            set(self.t2_now_edt5, 'unit', 'pixels', 'position', [900 503 70 23]);
            set(self.t2_now_edt5, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt5, 'Visible', 'off');
            set(self.t2_now_edt5, 'String', self.user_inputs.tab_2_now.midas_exogenous_lags);
            set(self.t2_now_edt5, 'CallBack', @self.cb_t2_now_edt5);

            % endogenous tightness label
            self.t2_now_txt12 = uicontrol('style', 'text');
            set(self.t2_now_txt12, 'unit', 'pixels', 'position', [520 471 200 25]);
            set(self.t2_now_txt12, 'String', ' ω₁:endo tightness');
            set(self.t2_now_txt12, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt12, 'FontName', 'Serif');
            set(self.t2_now_txt12, 'FontSize', 11);
            set(self.t2_now_txt12, 'FontWeight', 'bold');
            set(self.t2_now_txt12, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt12, 'Visible', 'off');  

            % endogenous tightness edit
            self.t2_now_edt7 = uicontrol('style','edit');
            set(self.t2_now_edt7, 'unit', 'pixels', 'position', [680 476 70 23]);
            set(self.t2_now_edt7, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt7, 'Visible', 'off');
            set(self.t2_now_edt7, 'String', self.user_inputs.tab_2_now.midas_omega1);
            set(self.t2_now_edt7, 'CallBack', @self.cb_t2_now_edt7);

            % endogenous lag decay label
            self.t2_now_txt13 = uicontrol('style', 'text');
            set(self.t2_now_txt13, 'unit', 'pixels', 'position', [520 446 200 25]);
            set(self.t2_now_txt13, 'String', ' ω₂:endo lag decay');
            set(self.t2_now_txt13, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt13, 'FontName', 'Serif');
            set(self.t2_now_txt13, 'FontSize', 11);
            set(self.t2_now_txt13, 'FontWeight', 'bold');
            set(self.t2_now_txt13, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt13, 'Visible', 'off');  

            % endogenous lag decay edit
            self.t2_now_edt8 = uicontrol('style','edit');
            set(self.t2_now_edt8, 'unit', 'pixels', 'position', [680 450 70 23]);
            set(self.t2_now_edt8, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt8, 'Visible', 'off');
            set(self.t2_now_edt8, 'String', self.user_inputs.tab_2_now.midas_omega2);
            set(self.t2_now_edt8, 'CallBack', @self.cb_t2_now_edt8);

            % exogenous tightness label
            self.t2_now_txt14 = uicontrol('style', 'text');
            set(self.t2_now_txt14, 'unit', 'pixels', 'position', [750 471 200 25]);
            set(self.t2_now_txt14, 'String', ' υ₁: exo tightness');
            set(self.t2_now_txt14, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt14, 'FontName', 'Serif');
            set(self.t2_now_txt14, 'FontSize', 11);
            set(self.t2_now_txt14, 'FontWeight', 'bold');
            set(self.t2_now_txt14, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt14, 'Visible', 'off');  

            % exogenous tightness edit
            self.t2_now_edt9 = uicontrol('style','edit');
            set(self.t2_now_edt9, 'unit', 'pixels', 'position', [900 476 70 23]);
            set(self.t2_now_edt9, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt9, 'Visible', 'off');
            set(self.t2_now_edt9, 'String', self.user_inputs.tab_2_now.midas_upsilon1);
            set(self.t2_now_edt9, 'CallBack', @self.cb_t2_now_edt9);

            % exogenous lag decay label
            self.t2_now_txt15 = uicontrol('style', 'text');
            set(self.t2_now_txt15, 'unit', 'pixels', 'position', [750 446 200 25]);
            set(self.t2_now_txt15, 'String', ' υ₂: exo lag decay');
            set(self.t2_now_txt15, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt15, 'FontName', 'Serif');
            set(self.t2_now_txt15, 'FontSize', 11);
            set(self.t2_now_txt15, 'FontWeight', 'bold');
            set(self.t2_now_txt15, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt15, 'Visible', 'off');  

            % exogenous lag decay edit
            self.t2_now_edt10 = uicontrol('style','edit');
            set(self.t2_now_edt10, 'unit', 'pixels', 'position', [900 450 70 23]);
            set(self.t2_now_edt10, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt10, 'Visible', 'off');
            set(self.t2_now_edt10, 'String', self.user_inputs.tab_2_now.midas_upsilon2);
            set(self.t2_now_edt10, 'CallBack', @self.cb_t2_now_edt10);

            % polynomial order label
            self.t2_now_txt10 = uicontrol('style', 'text');
            set(self.t2_now_txt10, 'unit', 'pixels', 'position', [520 424 200 20]);
            set(self.t2_now_txt10, 'String', ' polynomial order');
            set(self.t2_now_txt10, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt10, 'FontName', 'Serif');
            set(self.t2_now_txt10, 'FontSize', 11);
            set(self.t2_now_txt10, 'FontWeight', 'bold');
            set(self.t2_now_txt10, 'BackgroundColor', self.background_color); 
            set(self.t2_now_txt10, 'Visible', 'off');  

            % polynomial order edit
            self.t2_now_edt6 = uicontrol('style','edit');
            set(self.t2_now_edt6, 'unit', 'pixels', 'position', [680 425 70 23]);
            set(self.t2_now_edt6, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt6, 'Visible', 'off');
            set(self.t2_now_edt6, 'String', self.user_inputs.tab_2_now.midas_polynomial_order);
            set(self.t2_now_edt6, 'CallBack', @self.cb_t2_now_edt6);

            % mfbvar label
            self.t2_now_txt16 = uicontrol('style', 'text');
            set(self.t2_now_txt16, 'unit', 'pixels', 'position', [30 365 400 30]);
            set(self.t2_now_txt16, 'String', ' Mixed frequency Bayesian VAR');
            set(self.t2_now_txt16, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt16, 'FontName', 'Serif');
            set(self.t2_now_txt16, 'FontSize', 15);
            set(self.t2_now_txt16, 'FontWeight', 'bold');
            set(self.t2_now_txt16, 'FontAngle', 'italic');
            set(self.t2_now_txt16, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt16, 'Visible', 'off');
            
            % frame around mfbvar
            self.t2_now_frm3 = uicontrol('style','frame');
            set(self.t2_now_frm3, 'unit', 'pixels', 'position', [20 20 470 347]);
            set(self.t2_now_frm3, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_now_frm3, 'BackgroundColor', self.background_color);
            set(self.t2_now_frm3, 'Visible', 'off'); 

            % constant label
            self.t2_now_txt17 = uicontrol('style', 'text');
            set(self.t2_now_txt17, 'unit', 'pixels', 'position', [30 325 400 30]);
            set(self.t2_now_txt17, 'String', ' constant');
            set(self.t2_now_txt17, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt17, 'FontName', 'Serif');
            set(self.t2_now_txt17, 'FontSize', 11);
            set(self.t2_now_txt17, 'FontWeight', 'bold');
            set(self.t2_now_txt17, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt17, 'Visible', 'off');

            % constant checkbox
            self.t2_now_cbx1 = uicontrol('style', 'checkbox');
            set(self.t2_now_cbx1, 'unit', 'pixels', 'position', [460 335 20 20]);
            set(self.t2_now_cbx1, 'BackgroundColor', self.background_color);
            set(self.t2_now_cbx1, 'Visible', 'off');
            set(self.t2_now_cbx1, 'Value', self.user_inputs.tab_2_now.mfbvar_constant);
            set(self.t2_now_cbx1, 'CallBack', @self.cb_t2_now_cbx1);

            % linear trend label
            self.t2_now_txt18 = uicontrol('style', 'text');
            set(self.t2_now_txt18, 'unit', 'pixels', 'position', [30 298 400 30]);
            set(self.t2_now_txt18, 'String', ' linear trend');
            set(self.t2_now_txt18, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt18, 'FontName', 'Serif');
            set(self.t2_now_txt18, 'FontSize', 11);
            set(self.t2_now_txt18, 'FontWeight', 'bold');
            set(self.t2_now_txt18, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt18, 'Visible', 'off');             

            % linear trend checkbox
            self.t2_now_cbx2 = uicontrol('style', 'checkbox');
            set(self.t2_now_cbx2, 'unit', 'pixels', 'position', [460 308 20 20]);
            set(self.t2_now_cbx2, 'BackgroundColor', self.background_color);
            set(self.t2_now_cbx2, 'Visible', 'off');
            set(self.t2_now_cbx2, 'Value', self.user_inputs.tab_2_now.mfbvar_trend);
            set(self.t2_now_cbx2, 'CallBack', @self.cb_t2_now_cbx2);

            % quadratic trend label
            self.t2_now_txt19 = uicontrol('style', 'text');
            set(self.t2_now_txt19, 'unit', 'pixels', 'position', [30 271 400 30]);
            set(self.t2_now_txt19, 'String', ' quadratic trend');
            set(self.t2_now_txt19, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt19, 'FontName', 'Serif');
            set(self.t2_now_txt19, 'FontSize', 11);
            set(self.t2_now_txt19, 'FontWeight', 'bold');
            set(self.t2_now_txt19, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt19, 'Visible', 'off'); 

            % quadratic trend checkbox
            self.t2_now_cbx3 = uicontrol('style', 'checkbox');
            set(self.t2_now_cbx3, 'unit', 'pixels', 'position', [460 281 20 20]);
            set(self.t2_now_cbx3, 'BackgroundColor', self.background_color);
            set(self.t2_now_cbx3, 'Visible', 'off');
            set(self.t2_now_cbx3, 'Value', self.user_inputs.tab_2_now.mfbvar_quadratic_trend);
            set(self.t2_now_cbx3, 'CallBack', @self.cb_t2_now_cbx3);

            % decomposition label
            self.t2_now_txt20 = uicontrol('style', 'text');
            set(self.t2_now_txt20, 'unit', 'pixels', 'position', [30 244 400 30]);
            set(self.t2_now_txt20, 'String', ' decomposition');
            set(self.t2_now_txt20, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt20, 'FontName', 'Serif');
            set(self.t2_now_txt20, 'FontSize', 11);
            set(self.t2_now_txt20, 'FontWeight', 'bold');
            set(self.t2_now_txt20, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt20, 'Visible', 'off'); 

            % decomposition checkbox
            self.t2_now_cbx4 = uicontrol('style', 'checkbox');
            set(self.t2_now_cbx4, 'unit', 'pixels', 'position', [460 254 20 20]);
            set(self.t2_now_cbx4, 'BackgroundColor', self.background_color);
            set(self.t2_now_cbx4, 'Visible', 'off');
            set(self.t2_now_cbx4, 'Value', self.user_inputs.tab_2_now.mfbvar_decomposition);
            set(self.t2_now_cbx4, 'CallBack', @self.cb_t2_now_cbx4);

            % lags label
            self.t2_now_txt21 = uicontrol('style', 'text');
            set(self.t2_now_txt21, 'unit', 'pixels', 'position', [30 217 400 30]);
            set(self.t2_now_txt21, 'String', ' p:    lags');
            set(self.t2_now_txt21, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt21, 'FontName', 'Serif');
            set(self.t2_now_txt21, 'FontSize', 11);
            set(self.t2_now_txt21, 'FontWeight', 'bold');
            set(self.t2_now_txt21, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt21, 'Visible', 'off'); 

            % lags edit
            self.t2_now_edt11 = uicontrol('style','edit');
            set(self.t2_now_edt11, 'unit', 'pixels', 'position', [333 225 140 23]);
            set(self.t2_now_edt11, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt11, 'Visible', 'off');
            set(self.t2_now_edt11, 'String', self.user_inputs.tab_2_now.mfbvar_lags);
            set(self.t2_now_edt11, 'CallBack', @self.cb_t2_now_edt11);

            % ar coefficients label
            self.t2_now_txt22 = uicontrol('style', 'text');
            set(self.t2_now_txt22, 'unit', 'pixels', 'position', [30 190 400 30]);
            set(self.t2_now_txt22, 'String', ' δ:    AR coefficients');
            set(self.t2_now_txt22, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt22, 'FontName', 'Serif');
            set(self.t2_now_txt22, 'FontSize', 11);
            set(self.t2_now_txt22, 'FontWeight', 'bold');
            set(self.t2_now_txt22, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt22, 'Visible', 'off'); 

            % ar coefficients edit
            self.t2_now_edt12 = uicontrol('style','edit');
            set(self.t2_now_edt12, 'unit', 'pixels', 'position', [333 198 140 23]);
            set(self.t2_now_edt12, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt12, 'Visible', 'off');
            set(self.t2_now_edt12, 'String', self.user_inputs.tab_2_now.mfbvar_ar_coefficients);
            set(self.t2_now_edt12, 'CallBack', @self.cb_t2_now_edt12);

            % pi1 label
            self.t2_now_txt23 = uicontrol('style', 'text');
            set(self.t2_now_txt23, 'unit', 'pixels', 'position', [30 163 400 30]);
            set(self.t2_now_txt23, 'String', ' π₁:  overall tightness');
            set(self.t2_now_txt23, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt23, 'FontName', 'Serif');
            set(self.t2_now_txt23, 'FontSize', 11);
            set(self.t2_now_txt23, 'FontWeight', 'bold');
            set(self.t2_now_txt23, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt23, 'Visible', 'off');  

            % pi1 edit
            self.t2_now_edt13 = uicontrol('style','edit');
            set(self.t2_now_edt13, 'unit', 'pixels', 'position', [333 171 140 23]);
            set(self.t2_now_edt13, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt13, 'Visible', 'off');
            set(self.t2_now_edt13, 'String', self.user_inputs.tab_2_now.mfbvar_pi1);
            set(self.t2_now_edt13, 'CallBack', @self.cb_t2_now_edt13);

            % pi2 label
            self.t2_now_txt24 = uicontrol('style', 'text');
            set(self.t2_now_txt24, 'unit', 'pixels', 'position', [30 136 400 30]);
            set(self.t2_now_txt24, 'String', ' π₂:  cross-variable shrinkage');
            set(self.t2_now_txt24, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt24, 'FontName', 'Serif');
            set(self.t2_now_txt24, 'FontSize', 11);
            set(self.t2_now_txt24, 'FontWeight', 'bold');
            set(self.t2_now_txt24, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt24, 'Visible', 'off'); 

            % pi2 edit
            self.t2_now_edt14 = uicontrol('style','edit');
            set(self.t2_now_edt14, 'unit', 'pixels', 'position', [333 144 140 23]);
            set(self.t2_now_edt14, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt14, 'Visible', 'off');
            set(self.t2_now_edt14, 'String', self.user_inputs.tab_2_now.mfbvar_pi2);
            set(self.t2_now_edt14, 'CallBack', @self.cb_t2_now_edt14);

            % pi3 label
            self.t2_now_txt25 = uicontrol('style', 'text');
            set(self.t2_now_txt25, 'unit', 'pixels', 'position', [30 109 400 30]);
            set(self.t2_now_txt25, 'String', ' π₃:  lag decay');
            set(self.t2_now_txt25, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt25, 'FontName', 'Serif');
            set(self.t2_now_txt25, 'FontSize', 11);
            set(self.t2_now_txt25, 'FontWeight', 'bold');
            set(self.t2_now_txt25, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt25, 'Visible', 'off'); 

            % pi3 edit
            self.t2_now_edt15 = uicontrol('style','edit');
            set(self.t2_now_edt15, 'unit', 'pixels', 'position', [333 117 140 23]);
            set(self.t2_now_edt15, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt15, 'Visible', 'off');
            set(self.t2_now_edt15, 'String', self.user_inputs.tab_2_now.mfbvar_pi3);
            set(self.t2_now_edt15, 'CallBack', @self.cb_t2_now_edt15);

            % pi4 label
            self.t2_now_txt26 = uicontrol('style', 'text');
            set(self.t2_now_txt26, 'unit', 'pixels', 'position', [30 82 400 30]);
            set(self.t2_now_txt26, 'String', ' π₄:  exogenous slackness');
            set(self.t2_now_txt26, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt26, 'FontName', 'Serif');
            set(self.t2_now_txt26, 'FontSize', 11);
            set(self.t2_now_txt26, 'FontWeight', 'bold');
            set(self.t2_now_txt26, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt26, 'Visible', 'off'); 

            % pi4 edit
            self.t2_now_edt16 = uicontrol('style','edit');
            set(self.t2_now_edt16, 'unit', 'pixels', 'position', [333 90 140 23]);
            set(self.t2_now_edt16, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt16, 'Visible', 'off');
            set(self.t2_now_edt16, 'String', self.user_inputs.tab_2_now.mfbvar_pi4);
            set(self.t2_now_edt16, 'CallBack', @self.cb_t2_now_edt16);

            % decomposition file label
            self.t2_now_txt27 = uicontrol('style', 'text');
            set(self.t2_now_txt27, 'unit', 'pixels', 'position', [30 55 400 30]);
            set(self.t2_now_txt27, 'String', ' file: decomposition table');
            set(self.t2_now_txt27, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt27, 'FontName', 'Serif');
            set(self.t2_now_txt27, 'FontSize', 11);
            set(self.t2_now_txt27, 'FontWeight', 'bold');
            set(self.t2_now_txt27, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt27, 'Visible', 'off');  

            % decomposition file edit
            self.t2_now_edt17 = uicontrol('style','edit');
            set(self.t2_now_edt17, 'unit', 'pixels', 'position', [33 35 440 23]);
            set(self.t2_now_edt17, 'HorizontalAlignment', 'left');  
            set(self.t2_now_edt17, 'Visible', 'off');
            set(self.t2_now_edt17, 'String', self.user_inputs.tab_2_now.mfbvar_decomposition_file);
            set(self.t2_now_edt17, 'CallBack', @self.cb_t2_now_edt17);

            % dfm label
            self.t2_now_txt28 = uicontrol('style', 'text');
            set(self.t2_now_txt28, 'unit', 'pixels', 'position', [520 365 400 30]);
            set(self.t2_now_txt28, 'String', ' Bayesian dynamic factor model');
            set(self.t2_now_txt28, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt28, 'FontName', 'Serif');
            set(self.t2_now_txt28, 'FontSize', 15);
            set(self.t2_now_txt28, 'FontWeight', 'bold');
            set(self.t2_now_txt28, 'FontAngle', 'italic');
            set(self.t2_now_txt28, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt28, 'Visible', 'off');
            
            % frame around dfm
            self.t2_now_frm4 = uicontrol('style','frame');
            set(self.t2_now_frm4, 'unit', 'pixels', 'position', [510 20 470 347]);
            set(self.t2_now_frm4, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_now_frm4, 'BackgroundColor', self.background_color);
            set(self.t2_now_frm4, 'Visible', 'off'); 

            % factor label
            self.t2_now_txt29 = uicontrol('style', 'text');
            set(self.t2_now_txt29, 'unit', 'pixels', 'position', [520 325 400 30]);
            set(self.t2_now_txt29, 'String', ' m:   factors');
            set(self.t2_now_txt29, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt29, 'FontName', 'Serif');
            set(self.t2_now_txt29, 'FontSize', 11);
            set(self.t2_now_txt29, 'FontWeight', 'bold');
            set(self.t2_now_txt29, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt29, 'Visible', 'off');

            % factor edit
            self.t2_now_edt18 = uicontrol('style','edit');
            set(self.t2_now_edt18, 'unit', 'pixels', 'position', [825 333 140 23]);
            set(self.t2_now_edt18, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt18, 'Visible', 'off');
            set(self.t2_now_edt18, 'String', self.user_inputs.tab_2_now.dfm_factors);
            set(self.t2_now_edt18, 'CallBack', @self.cb_t2_now_edt18);

            % loadings lag label
            self.t2_now_txt30 = uicontrol('style', 'text');
            set(self.t2_now_txt30, 'unit', 'pixels', 'position', [520 298 400 30]);
            set(self.t2_now_txt30, 'String', ' q:    loadings lags');
            set(self.t2_now_txt30, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt30, 'FontName', 'Serif');
            set(self.t2_now_txt30, 'FontSize', 11);
            set(self.t2_now_txt30, 'FontWeight', 'bold');
            set(self.t2_now_txt30, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt30, 'Visible', 'off');

            % loadings lag edit
            self.t2_now_edt19 = uicontrol('style','edit');
            set(self.t2_now_edt19, 'unit', 'pixels', 'position', [825 306 140 23]);
            set(self.t2_now_edt19, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt19, 'Visible', 'off');
            set(self.t2_now_edt19, 'String', self.user_inputs.tab_2_now.dfm_loadings_lags);
            set(self.t2_now_edt19, 'CallBack', @self.cb_t2_now_edt19);

            % factor lag label
            self.t2_now_txt31 = uicontrol('style', 'text');
            set(self.t2_now_txt31, 'unit', 'pixels', 'position', [520 271 400 30]);
            set(self.t2_now_txt31, 'String', ' p:    factor lags');
            set(self.t2_now_txt31, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt31, 'FontName', 'Serif');
            set(self.t2_now_txt31, 'FontSize', 11);
            set(self.t2_now_txt31, 'FontWeight', 'bold');
            set(self.t2_now_txt31, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt31, 'Visible', 'off');

            % factor lag edit
            self.t2_now_edt20 = uicontrol('style','edit');
            set(self.t2_now_edt20, 'unit', 'pixels', 'position', [825 279 140 23]);
            set(self.t2_now_edt20, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt20, 'Visible', 'off');
            set(self.t2_now_edt20, 'String', self.user_inputs.tab_2_now.dfm_factor_lags);
            set(self.t2_now_edt20, 'CallBack', @self.cb_t2_now_edt20);

            % residual lag label
            self.t2_now_txt32 = uicontrol('style', 'text');
            set(self.t2_now_txt32, 'unit', 'pixels', 'position', [520 244 400 30]);
            set(self.t2_now_txt32, 'String', ' r:    residual lags');
            set(self.t2_now_txt32, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt32, 'FontName', 'Serif');
            set(self.t2_now_txt32, 'FontSize', 11);
            set(self.t2_now_txt32, 'FontWeight', 'bold');
            set(self.t2_now_txt32, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt32, 'Visible', 'off');

            % residual lag edit
            self.t2_now_edt21 = uicontrol('style','edit');
            set(self.t2_now_edt21, 'unit', 'pixels', 'position', [825 252 140 23]);
            set(self.t2_now_edt21, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt21, 'Visible', 'off');
            set(self.t2_now_edt21, 'String', self.user_inputs.tab_2_now.dfm_residual_lags);
            set(self.t2_now_edt21, 'CallBack', @self.cb_t2_now_edt21);

            % residual variance label
            self.t2_now_txt33 = uicontrol('style', 'text');
            set(self.t2_now_txt33, 'unit', 'pixels', 'position', [520 217 400 30]);
            set(self.t2_now_txt33, 'String', ' σ:    residual variance');
            set(self.t2_now_txt33, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt33, 'FontName', 'Serif');
            set(self.t2_now_txt33, 'FontSize', 11);
            set(self.t2_now_txt33, 'FontWeight', 'bold');
            set(self.t2_now_txt33, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt33, 'Visible', 'off');

            % residual variance edit
            self.t2_now_edt22 = uicontrol('style','edit');
            set(self.t2_now_edt22, 'unit', 'pixels', 'position', [825 225 140 23]);
            set(self.t2_now_edt22, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt22, 'Visible', 'off');
            set(self.t2_now_edt22, 'String', self.user_inputs.tab_2_now.dfm_sigma);
            set(self.t2_now_edt22, 'CallBack', @self.cb_t2_now_edt22);

            % factor variance label
            self.t2_now_txt34 = uicontrol('style', 'text');
            set(self.t2_now_txt34, 'unit', 'pixels', 'position', [520 190 400 30]);
            set(self.t2_now_txt34, 'String', ' ω:    factor variance');
            set(self.t2_now_txt34, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt34, 'FontName', 'Serif');
            set(self.t2_now_txt34, 'FontSize', 11);
            set(self.t2_now_txt34, 'FontWeight', 'bold');
            set(self.t2_now_txt34, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt34, 'Visible', 'off');

            % factor variance edit
            self.t2_now_edt23 = uicontrol('style','edit');
            set(self.t2_now_edt23, 'unit', 'pixels', 'position', [825 198 140 23]);
            set(self.t2_now_edt23, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt23, 'Visible', 'off');
            set(self.t2_now_edt23, 'String', self.user_inputs.tab_2_now.dfm_omega);
            set(self.t2_now_edt23, 'CallBack', @self.cb_t2_now_edt23);

            % loadings tightness label
            self.t2_now_txt35 = uicontrol('style', 'text');
            set(self.t2_now_txt35, 'unit', 'pixels', 'position', [520 163 400 30]);
            set(self.t2_now_txt35, 'String', ' δ₁:   loadings tightness');
            set(self.t2_now_txt35, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt35, 'FontName', 'Serif');
            set(self.t2_now_txt35, 'FontSize', 11);
            set(self.t2_now_txt35, 'FontWeight', 'bold');
            set(self.t2_now_txt35, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt35, 'Visible', 'off');

            % loadings tightness edit
            self.t2_now_edt24 = uicontrol('style','edit');
            set(self.t2_now_edt24, 'unit', 'pixels', 'position', [825 171 140 23]);
            set(self.t2_now_edt24, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt24, 'Visible', 'off');
            set(self.t2_now_edt24, 'String', self.user_inputs.tab_2_now.dfm_delta1);
            set(self.t2_now_edt24, 'CallBack', @self.cb_t2_now_edt24);

            % pi1 label
            self.t2_now_txt36 = uicontrol('style', 'text');
            set(self.t2_now_txt36, 'unit', 'pixels', 'position', [520 136 400 30]);
            set(self.t2_now_txt36, 'String', ' π₁:   factor tightness');
            set(self.t2_now_txt36, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt36, 'FontName', 'Serif');
            set(self.t2_now_txt36, 'FontSize', 11);
            set(self.t2_now_txt36, 'FontWeight', 'bold');
            set(self.t2_now_txt36, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt36, 'Visible', 'off');

            % pi1 edit
            self.t2_now_edt25 = uicontrol('style','edit');
            set(self.t2_now_edt25, 'unit', 'pixels', 'position', [825 144 140 23]);
            set(self.t2_now_edt25, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt25, 'Visible', 'off');
            set(self.t2_now_edt25, 'String', self.user_inputs.tab_2_now.dfm_pi1);
            set(self.t2_now_edt25, 'CallBack', @self.cb_t2_now_edt25);

            % pi2 label
            self.t2_now_txt37 = uicontrol('style', 'text');
            set(self.t2_now_txt37, 'unit', 'pixels', 'position', [520 109 400 30]);
            set(self.t2_now_txt37, 'String', ' π₂:   cross-variable shrinkage');
            set(self.t2_now_txt37, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt37, 'FontName', 'Serif');
            set(self.t2_now_txt37, 'FontSize', 11);
            set(self.t2_now_txt37, 'FontWeight', 'bold');
            set(self.t2_now_txt37, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt37, 'Visible', 'off');

            % pi2 edit
            self.t2_now_edt26 = uicontrol('style','edit');
            set(self.t2_now_edt26, 'unit', 'pixels', 'position', [825 117 140 23]);
            set(self.t2_now_edt26, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt26, 'Visible', 'off');
            set(self.t2_now_edt26, 'String', self.user_inputs.tab_2_now.dfm_pi2);
            set(self.t2_now_edt26, 'CallBack', @self.cb_t2_now_edt26);

            % pi3 label
            self.t2_now_txt38 = uicontrol('style', 'text');
            set(self.t2_now_txt38, 'unit', 'pixels', 'position', [520 82 400 30]);
            set(self.t2_now_txt38, 'String', ' π₃:   lag decay');
            set(self.t2_now_txt38, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt38, 'FontName', 'Serif');
            set(self.t2_now_txt38, 'FontSize', 11);
            set(self.t2_now_txt38, 'FontWeight', 'bold');
            set(self.t2_now_txt38, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt38, 'Visible', 'off');

            % pi3 edit
            self.t2_now_edt27 = uicontrol('style','edit');
            set(self.t2_now_edt27, 'unit', 'pixels', 'position', [825 90 140 23]);
            set(self.t2_now_edt27, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt27, 'Visible', 'off');
            set(self.t2_now_edt27, 'String', self.user_inputs.tab_2_now.dfm_pi3);
            set(self.t2_now_edt27, 'CallBack', @self.cb_t2_now_edt27);

            % omega1 label
            self.t2_now_txt39 = uicontrol('style', 'text');
            set(self.t2_now_txt39, 'unit', 'pixels', 'position', [520 55 400 30]);
            set(self.t2_now_txt39, 'String', ' ω₁:   residual tightness');
            set(self.t2_now_txt39, 'HorizontalAlignment', 'left');
            set(self.t2_now_txt39, 'FontName', 'Serif');
            set(self.t2_now_txt39, 'FontSize', 11);
            set(self.t2_now_txt39, 'FontWeight', 'bold');
            set(self.t2_now_txt39, 'BackgroundColor', self.background_color);  
            set(self.t2_now_txt39, 'Visible', 'off');

            % omega1 edit
            self.t2_now_edt28 = uicontrol('style','edit');
            set(self.t2_now_edt28, 'unit', 'pixels', 'position', [825 60 140 23]);
            set(self.t2_now_edt28, 'HorizontalAlignment', 'center');  
            set(self.t2_now_edt28, 'Visible', 'off');
            set(self.t2_now_edt28, 'String', self.user_inputs.tab_2_now.dfm_omega1);
            set(self.t2_now_edt28, 'CallBack', @self.cb_t2_now_edt28);
        end


        function hide_tab_2_now(self)

            % hide all controls            
            set(self.t2_now_txt1, 'Visible', 'off'); 
            set(self.t2_now_txt2, 'Visible', 'off');
            set(self.t2_now_txt3, 'Visible', 'off');
            set(self.t2_now_txt4, 'Visible', 'off');
            set(self.t2_now_txt5, 'Visible', 'off');
            set(self.t2_now_txt6, 'Visible', 'off');
            set(self.t2_now_txt7, 'Visible', 'off');
            set(self.t2_now_txt8, 'Visible', 'off');
            set(self.t2_now_txt9, 'Visible', 'off');
            set(self.t2_now_txt10, 'Visible', 'off');
            set(self.t2_now_txt11, 'Visible', 'off');
            set(self.t2_now_txt12, 'Visible', 'off');
            set(self.t2_now_txt13, 'Visible', 'off');
            set(self.t2_now_txt14, 'Visible', 'off');
            set(self.t2_now_txt15, 'Visible', 'off');
            set(self.t2_now_txt16, 'Visible', 'off');
            set(self.t2_now_txt17, 'Visible', 'off');
            set(self.t2_now_txt18, 'Visible', 'off');
            set(self.t2_now_txt19, 'Visible', 'off');
            set(self.t2_now_txt20, 'Visible', 'off');
            set(self.t2_now_txt21, 'Visible', 'off');
            set(self.t2_now_txt22, 'Visible', 'off');
            set(self.t2_now_txt23, 'Visible', 'off');
            set(self.t2_now_txt24, 'Visible', 'off');
            set(self.t2_now_txt25, 'Visible', 'off');
            set(self.t2_now_txt26, 'Visible', 'off');
            set(self.t2_now_txt27, 'Visible', 'off');
            set(self.t2_now_txt28, 'Visible', 'off');
            set(self.t2_now_txt29, 'Visible', 'off');
            set(self.t2_now_txt30, 'Visible', 'off');
            set(self.t2_now_txt31, 'Visible', 'off');
            set(self.t2_now_txt32, 'Visible', 'off');
            set(self.t2_now_txt33, 'Visible', 'off');
            set(self.t2_now_txt34, 'Visible', 'off');
            set(self.t2_now_txt35, 'Visible', 'off');
            set(self.t2_now_txt36, 'Visible', 'off');
            set(self.t2_now_txt37, 'Visible', 'off');
            set(self.t2_now_txt38, 'Visible', 'off');
            set(self.t2_now_txt39, 'Visible', 'off');
            set(self.t2_now_frm1, 'Visible', 'off');
            set(self.t2_now_frm2, 'Visible', 'off');
            set(self.t2_now_frm3, 'Visible', 'off');
            set(self.t2_now_frm4, 'Visible', 'off');
            set(self.t2_now_bgr1, 'Visible', 'off');
            set(self.t2_now_rdb1, 'Visible', 'off');
            set(self.t2_now_rdb2, 'Visible', 'off');
            set(self.t2_now_rdb3, 'Visible', 'off');
            set(self.t2_now_edt1, 'Visible', 'off');  
            set(self.t2_now_edt2, 'Visible', 'off');  
            set(self.t2_now_edt3, 'Visible', 'off');  
            set(self.t2_now_edt4, 'Visible', 'off');  
            set(self.t2_now_edt5, 'Visible', 'off');  
            set(self.t2_now_edt6, 'Visible', 'off');  
            set(self.t2_now_edt7, 'Visible', 'off');  
            set(self.t2_now_edt8, 'Visible', 'off');
            set(self.t2_now_edt9, 'Visible', 'off');
            set(self.t2_now_edt10, 'Visible', 'off');
            set(self.t2_now_edt11, 'Visible', 'off');
            set(self.t2_now_edt12, 'Visible', 'off');
            set(self.t2_now_edt13, 'Visible', 'off');
            set(self.t2_now_edt14, 'Visible', 'off');
            set(self.t2_now_edt15, 'Visible', 'off');
            set(self.t2_now_edt16, 'Visible', 'off');
            set(self.t2_now_edt17, 'Visible', 'off');
            set(self.t2_now_edt18, 'Visible', 'off');
            set(self.t2_now_edt19, 'Visible', 'off');
            set(self.t2_now_edt20, 'Visible', 'off');
            set(self.t2_now_edt21, 'Visible', 'off');
            set(self.t2_now_edt22, 'Visible', 'off');
            set(self.t2_now_edt23, 'Visible', 'off');
            set(self.t2_now_edt24, 'Visible', 'off');
            set(self.t2_now_edt25, 'Visible', 'off');
            set(self.t2_now_edt26, 'Visible', 'off');
            set(self.t2_now_edt27, 'Visible', 'off');
            set(self.t2_now_edt28, 'Visible', 'off');
            set(self.t2_now_mnu1, 'Visible', 'off');
            set(self.t2_now_cbx1, 'Visible', 'off');
            set(self.t2_now_cbx2, 'Visible', 'off');
            set(self.t2_now_cbx3, 'Visible', 'off');
            set(self.t2_now_cbx4, 'Visible', 'off');

            % update tab color
            set(self.tab_pbt2, 'BackgroundColor', self.backtabs_color);

        end


        function show_tab_2_now(self)

            % show all controls
            set(self.t2_now_txt1, 'Visible', 'on'); 
            set(self.t2_now_txt2, 'Visible', 'on'); 
            set(self.t2_now_txt3, 'Visible', 'on'); 
            set(self.t2_now_txt4, 'Visible', 'on'); 
            set(self.t2_now_txt5, 'Visible', 'on'); 
            set(self.t2_now_txt6, 'Visible', 'on'); 
            set(self.t2_now_txt7, 'Visible', 'on'); 
            set(self.t2_now_txt8, 'Visible', 'on'); 
            set(self.t2_now_txt9, 'Visible', 'on'); 
            set(self.t2_now_txt10, 'Visible', 'on'); 
            set(self.t2_now_txt11, 'Visible', 'on'); 
            set(self.t2_now_txt12, 'Visible', 'on'); 
            set(self.t2_now_txt13, 'Visible', 'on'); 
            set(self.t2_now_txt14, 'Visible', 'on'); 
            set(self.t2_now_txt15, 'Visible', 'on'); 
            set(self.t2_now_txt16, 'Visible', 'on'); 
            set(self.t2_now_txt17, 'Visible', 'on'); 
            set(self.t2_now_txt18, 'Visible', 'on'); 
            set(self.t2_now_txt19, 'Visible', 'on'); 
            set(self.t2_now_txt20, 'Visible', 'on'); 
            set(self.t2_now_txt21, 'Visible', 'on'); 
            set(self.t2_now_txt22, 'Visible', 'on'); 
            set(self.t2_now_txt23, 'Visible', 'on'); 
            set(self.t2_now_txt24, 'Visible', 'on'); 
            set(self.t2_now_txt25, 'Visible', 'on'); 
            set(self.t2_now_txt26, 'Visible', 'on'); 
            set(self.t2_now_txt27, 'Visible', 'on'); 
            set(self.t2_now_txt28, 'Visible', 'on'); 
            set(self.t2_now_txt29, 'Visible', 'on'); 
            set(self.t2_now_txt30, 'Visible', 'on'); 
            set(self.t2_now_txt31, 'Visible', 'on'); 
            set(self.t2_now_txt32, 'Visible', 'on'); 
            set(self.t2_now_txt33, 'Visible', 'on'); 
            set(self.t2_now_txt34, 'Visible', 'on'); 
            set(self.t2_now_txt35, 'Visible', 'on'); 
            set(self.t2_now_txt36, 'Visible', 'on'); 
            set(self.t2_now_txt37, 'Visible', 'on'); 
            set(self.t2_now_txt38, 'Visible', 'on'); 
            set(self.t2_now_txt39, 'Visible', 'on'); 
            set(self.t2_now_frm1, 'Visible', 'on');
            set(self.t2_now_frm2, 'Visible', 'on');
            set(self.t2_now_frm3, 'Visible', 'on');
            set(self.t2_now_frm4, 'Visible', 'on');
            set(self.t2_now_bgr1, 'Visible', 'on');
            set(self.t2_now_rdb1, 'Visible', 'on');
            set(self.t2_now_rdb2, 'Visible', 'on');
            set(self.t2_now_rdb3, 'Visible', 'on');
            set(self.t2_now_edt1, 'Visible', 'on');    
            set(self.t2_now_edt2, 'Visible', 'on');  
            set(self.t2_now_edt3, 'Visible', 'on'); 
            set(self.t2_now_edt4, 'Visible', 'on'); 
            set(self.t2_now_edt5, 'Visible', 'on'); 
            set(self.t2_now_edt6, 'Visible', 'on'); 
            set(self.t2_now_edt7, 'Visible', 'on'); 
            set(self.t2_now_edt8, 'Visible', 'on'); 
            set(self.t2_now_edt9, 'Visible', 'on'); 
            set(self.t2_now_edt10, 'Visible', 'on'); 
            set(self.t2_now_edt11, 'Visible', 'on'); 
            set(self.t2_now_edt12, 'Visible', 'on'); 
            set(self.t2_now_edt13, 'Visible', 'on'); 
            set(self.t2_now_edt14, 'Visible', 'on'); 
            set(self.t2_now_edt15, 'Visible', 'on'); 
            set(self.t2_now_edt16, 'Visible', 'on'); 
            set(self.t2_now_edt17, 'Visible', 'on'); 
            set(self.t2_now_edt18, 'Visible', 'on'); 
            set(self.t2_now_edt19, 'Visible', 'on'); 
            set(self.t2_now_edt20, 'Visible', 'on'); 
            set(self.t2_now_edt21, 'Visible', 'on'); 
            set(self.t2_now_edt22, 'Visible', 'on'); 
            set(self.t2_now_edt23, 'Visible', 'on'); 
            set(self.t2_now_edt24, 'Visible', 'on'); 
            set(self.t2_now_edt25, 'Visible', 'on'); 
            set(self.t2_now_edt26, 'Visible', 'on'); 
            set(self.t2_now_edt27, 'Visible', 'on'); 
            set(self.t2_now_edt28, 'Visible', 'on'); 
            set(self.t2_now_mnu1, 'Visible', 'on');
            set(self.t2_now_cbx1, 'Visible', 'on');
            set(self.t2_now_cbx2, 'Visible', 'on');
            set(self.t2_now_cbx3, 'Visible', 'on');
            set(self.t2_now_cbx4, 'Visible', 'on');
        end
        
        function cb_t2_now_bgr1(self, hObject, callbackdata)
           if get(self.t2_now_bgr1, 'SelectedObject') == self.t2_now_rdb1
               self.user_inputs.tab_2_now.model = 1;
           elseif get(self.t2_now_bgr1, 'SelectedObject') == self.t2_now_rdb2
               self.user_inputs.tab_2_now.model = 2;
           elseif get(self.t2_now_bgr1, 'SelectedObject') == self.t2_now_rdb3
               self.user_inputs.tab_2_now.model = 3;               
           end
        end

        function cb_t2_now_edt1(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.iterations = get(self.t2_now_edt1, 'String');
        end 

        function cb_t2_now_edt2(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.burnin = get(self.t2_now_edt2, 'String');
        end 

        function cb_t2_now_edt3(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.model_credibility = get(self.t2_now_edt3, 'String');
        end 

        function cb_t2_now_edt4(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_endogenous_lags = get(self.t2_now_edt4, 'String');
        end 

        function cb_t2_now_edt5(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_exogenous_lags = get(self.t2_now_edt5, 'String');
        end 

        function cb_t2_now_edt6(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_polynomial_order = get(self.t2_now_edt6, 'String');
        end 

        function cb_t2_now_mnu1(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_model = get(self.t2_now_mnu1, 'Value');
            % self.user_inputs.tab_2_now.midas_prior_type = get(self.t2_now_mnu1, 'Value');
        end

        function cb_t2_now_edt7(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_omega1 = get(self.t2_now_edt7, 'String');
        end 

        function cb_t2_now_edt8(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_omega2 = get(self.t2_now_edt8, 'String');
        end 

        function cb_t2_now_edt9(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_upsilon1 = get(self.t2_now_edt9, 'String');
        end 

        function cb_t2_now_edt10(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.midas_upsilon2 = get(self.t2_now_edt10, 'String');
        end 

        function cb_t2_now_cbx1(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_constant = logical(get(self.t2_now_cbx1, 'Value'));
        end 

        function cb_t2_now_cbx2(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_trend = logical(get(self.t2_now_cbx2, 'Value'));
        end 

        function cb_t2_now_cbx3(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_quadratic_trend = logical(get(self.t2_now_cbx3, 'Value'));
        end 

        function cb_t2_now_cbx4(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_decomposition = logical(get(self.t2_now_cbx4, 'Value'));
        end 

        function cb_t2_now_edt11(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_lags = get(self.t2_now_edt11, 'String');
        end 

        function cb_t2_now_edt12(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_ar_coefficients = get(self.t2_now_edt12, 'String');
        end 

        function cb_t2_now_edt13(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_pi1 = get(self.t2_now_edt13, 'String');
        end 

        function cb_t2_now_edt14(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_pi2 = get(self.t2_now_edt14, 'String');
        end 

        function cb_t2_now_edt15(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_pi3 = get(self.t2_now_edt15, 'String');
        end         

        function cb_t2_now_edt16(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_pi4 = get(self.t2_now_edt16, 'String');
        end 

        function cb_t2_now_edt17(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.mfbvar_decomposition_file = get(self.t2_now_edt17, 'String');
        end 

        function cb_t2_now_edt18(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_factors = get(self.t2_now_edt18, 'String');
        end 

        function cb_t2_now_edt19(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_loadings_lags = get(self.t2_now_edt19, 'String');
        end 

        function cb_t2_now_edt20(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_factor_lags = get(self.t2_now_edt20, 'String');
        end 

        function cb_t2_now_edt21(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_residual_lags = get(self.t2_now_edt21, 'String');
        end 

        function cb_t2_now_edt22(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_sigma = get(self.t2_now_edt22, 'String');
        end 

        function cb_t2_now_edt23(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_omega = get(self.t2_now_edt23, 'String');
        end 

        function cb_t2_now_edt24(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_delta1 = get(self.t2_now_edt24, 'String');
        end 

        function cb_t2_now_edt25(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_pi1 = get(self.t2_now_edt25, 'String');
        end 

        function cb_t2_now_edt26(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_pi2 = get(self.t2_now_edt26, 'String');
        end 

        function cb_t2_now_edt27(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_pi3 = get(self.t2_now_edt27, 'String');
        end 

        function cb_t2_now_edt28(self, hObject, callbackdata)
            self.user_inputs.tab_2_now.dfm_omega1 = get(self.t2_now_edt28, 'String');
        end 

        
    end  
     
end