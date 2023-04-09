%---------------------------------------------------
% Graphical User Interface
%--------------------------------------------------- 


% clear workspace and console
clear all;
close all;
warning off;
clc;

% add all toolbox files to workspace
addpath(genpath(fileparts(which(mfilename))));

% create and run graphical user interface
gui = GraphicalUserInterface();

% stop execution until user presses Run button (or closes GUI manually)
waitfor(gui, 'press_run');

% if user closed window manually without validating interface, terminate application
if gui.user_interrupt
    return;
end

% recover user inputs
user_inputs = gui.user_inputs;


%---------------------------------------------------
% Main code
%--------------------------------------------------- 


model = user_inputs.tab_1.model;


% if model is linear regression, run main code for linear regression, and return model
if model == 1
    lr = linear_regression_main_code(user_inputs);
end



