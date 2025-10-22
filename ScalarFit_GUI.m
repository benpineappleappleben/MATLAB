function fig = ScalarFit_GUI()
% SCALARFIT_GUI  Opens a GUI wrapper for the ScalarFit fitting utility.
%
% DESCRIPTION -------------------------------------------------------------
% ScalarFit_GUI provides a user-friendly graphical interface for the
% ScalarFit algorithm. To launch the GUI, run this function without input
% arguments or press the Run button in the MATLAB Editor. For information
% about the fitting algorithm and the name-value pairs, refer to the help
% article for ScalarFit.
%
% USING THE GUI -----------------------------------------------------------
% After opening the GUI, first choose a functional form from the dropdown
% list. If your desired function is not listed, select "Custom (Define
% Function)" and enter it in the Custom box. Use the format @(c,x) f(c,x),
% where c is a vector of parameters and x is the independent variable. Leave
% the box empty if you do not need a custom function.
%
% For example, to fit a model of the form c(1) + c(2)*x^2, type:
%   @(c,x) c(1) + c(2)*x.^2
% Since this function has two parameters, the initial guess vector must
% also contain two elements. MATLAB performs vectorized operations during
% fitting, so ensure all expressions use elementwise operators (e.g., .^,
% .*, ./). Otherwise, errors may occur.
%
% Next, specify the path to the data file to be fitted. The file must be a
% two-column matrix, with the independent variable in the first column and
% the dependent variable in the second.
%
% Finally, provide an initial guess vector. This guess serves as the
% starting point for the fit. It need not be exact, but better guesses
% improve both speed and accuracy.
%
% For details on available name-value arguments and predefined functional
% forms, see the documentation for ScalarFit and the file:
%   "ScalarFit Functional Forms v1.0.1.pdf"
%
% -------------------------------------------------------------------------

% /////////////////////////
% Benjamin Van Schaick
% Version 1.0.0
% Date Modified: 11/15/2024
% /////////////////////////

% Figure Initialization

fig = figure('Color', [1 1 1], 'Position', [400 25 510 660], 'Units', 'pixels');

% Title Interface

title_text = uicontrol('Parent', fig, 'Style', 'text', 'String', 'ScalarFit', 'ForegroundColor', [1 1 1], 'FontSize', 45, 'BackgroundColor', '#222498', 'HorizontalAlignment', 'center', 'Position', [20 565 350 75]);
fit_button = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'FIT', 'Callback', @buttonpress, 'BackgroundColor', '#222498', 'FontSize', 12, 'ForegroundColor', [1 1 1], 'Position', [390 570 100 65]);
form_text = uicontrol('Parent', fig, 'Style', 'text', 'String', 'Form: ', 'HorizontalAlignment', 'right', 'FontSize', 10, 'BackgroundColor', [1 1 1], 'Position', [20 519 50 30]);
form_dropdown = uicontrol('Parent', fig, 'Style', 'popupmenu', 'String', {'Linear','Quadratic','Cubic','Polynomial','Rational','Exponential','Logarithmic','Logistic','Gaussian','Sinusoid','Series (n-Step)','Custom (Define Function)'}, 'Position', [80 520 170 30]);
custom_text = uicontrol('Parent', fig, 'Style', 'text', 'String', 'Custom: ', 'FontSize', 10, 'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1], 'Position', [270 519 60 30]);
custom_fcn_input = uicontrol('Parent', fig, 'Style', 'edit', 'HorizontalAlignment', 'left', 'Position', [335 525 135 30]);
data_text = uicontrol('Parent', fig, 'Style', 'text', 'String', 'Data Path: ', 'BackgroundColor', [1 1 1], 'HorizontalAlignment', 'right', 'Position', [40 485 50 30]);
data_matrix_path = uicontrol('Parent', fig, 'Style', 'edit', 'HorizontalAlignment', 'left', 'Position', [100 485 135 30]);
guess_text = uicontrol('Parent', fig, 'Style', 'text', 'String', 'Guess: ', 'FontSize', 10, 'BackgroundColor', [1 1 1], 'HorizontalAlignment', 'right', 'Position', [275 480 50 30]);
guess_input = uicontrol('Parent', fig, 'Style', 'edit', 'HorizontalAlignment', 'left', 'Position', [335 485 135 30]);
varargin_text = uicontrol('Parent', fig, 'Style', 'text', 'String', 'Name-Value Pairs: ', 'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1], 'Position', [30 445 60 30]);
varargin_input = uicontrol('Parent', fig, 'Style', 'edit', 'HorizontalAlignment', 'left', 'Position', [100 445 370 30]);
ax = axes('Parent', fig, 'Position', [0.078431372549 0.21969696969 0.843137254902 0.409090909]);
grid(ax,'on')
title(ax,'Data and Fit')
result_textbox = uicontrol('Parent', fig, 'Style', 'text', 'String', '', 'FontSize', 10, 'BackgroundColor', [0.8 0.8 0.8], 'Position', [40 75 430 40]);
error_text = uicontrol('Parent', fig, 'Style', 'text', 'String', 'Error: ', 'FontSize', 12, 'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1], 'Position', [40 15 50 35]);
error_result_box = uicontrol('Parent', fig, 'Style', 'text', 'String', ' ', 'Position', [100 20 135 35]);
r_squared_box = uicontrol('Parent', fig, 'Style', 'text', 'String', 'R² = ', 'FontSize', 12, 'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1], 'Position', [275 15 50 35]);
r_squared_result_box = uicontrol('Parent', fig, 'Style', 'text', 'String', ' ', 'Position', [335 20 135 35]);

    function buttonpress(src,event)
        delete(ax.Children) % Clears the old plot (if one exists)
        forms = get(form_dropdown,'String'); % Gives a cell array containing all form names
        form = get(form_dropdown,'Value'); % Gives the index selected in the dropdown
        form = forms{form}; % Finds the name of the chosen form given the index and list of all names
        if isequal(form,'Custom (Define Function)') % If custom is chosen, we will rename to 'Custom' to simplify
            form = 'Custom';
        end

        data = get(data_matrix_path,'String'); % 
        data = readmatrix(data);
        guess = get(guess_input,'String');
        guess = str2num(guess);
        v = get(varargin_input,'String');

        if ~isequal(v,[])
            v = textscan(v,'%s','Delimiter',',');
            v = v{1};
        else
            v = {};
        end

        if isequal(form,'Custom')
            v = cat(1,v,{'SpecFcn';str2func(get(custom_fcn_input,'String'))});
        end

        v = cat(1,v,{'GUI';'true'});

        c = ScalarFit(form,data,guess,v);
        set(error_result_box,'String',num2str(c(length(c))))
        TSS = sum((mean(data(:,2)) - data(:,2)).^2);
        set(r_squared_result_box,'String',string(vpa(1-(c(length(c))/TSS),10)))

        output_parameters = '';
        for j = 1:1:length(c)-1
            output_parameters = cat(2,output_parameters,['c(',num2str(j),')=',num2str(c(j)),'   ']);
        end
        set(result_textbox,'String',output_parameters)
    end

end