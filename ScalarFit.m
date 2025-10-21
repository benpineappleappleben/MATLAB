function p = ScalarFit(form,data,guess,varargin)
% DESCRIPTION--------------------------------------------------------------
% 
% ScalarFit is a function which is capable of fitting a single independent 
% variable (x) to a single dependent variable (y) using any specified
% functional form and the least-squares regression method.
% 
% ARGUMENTS----------------------------------------------------------------
%
% The first argument of ScalarFit is the functional form, which must be entered
% as a character string. This argument chooses the function that will be
% used to fit your data. A list of all functional forms can be found at the
% end of this help article as well as in the appropriate PDF.
% 
% The second argument of ScalarFit is a matrix containing the data to be
% fitted. This data matrix must be an n-by-2 matrix with the independent
% variable in the first column and the dependent variable in the second
% column. The data matrix can also be inputted as a 2-by-n matrix, in which
% case the matrix will be transformed. In this case, the first row will be
% the independent variable data and the second row will be the dependent
% variable data. A 2-by-2 matrix will be left unmodified.
% 
% The third argument of ScalarFit is a row vector containing an initial guess
% for each parameter in the functional form. The meaning of the guess
% vector depends on the functional form chosen; refer to the end of this
% help article or to "ScalarFit Functional Forms v1.0.0.pdf" to see what the
% proper format and interpretation of your guess vector is. If your guess
% vector is entered as a column vector, it will be transformed into a row
% vector.
%
% VARARGIN ARGUMENTS-------------------------------------------------------
% 
% All varargin inputs should be entered in a name-value format as follows:
% ScalarFit(form,data,guess,'Name',value,'Name',value,...)
%
% Below is a list of the supported name-value pairs along with what they
% change in the ScalarFit pipeline. ScalarFit does not require a certain number of
% varargin arguments to operate unless you are using a form that requires a
% name-value pair. For example, the 'Rational' functional form requires a
% name-value pair to operate.
%
% In alphabetical order, the supported name-value pairs are:
%
% (1) 'Bounds' → A row vector containing a repeat of three values. The
% first index in this repeat is the index of the guess vector that you'd
% like to impose boundaries on. The second index is the lower bound, and
% the third index is the upper bound. You can impose boundaries on as many
% variables as you'd like, but the same trinumeric pattern must be kept.
% This means that the length of the vector will be a multiple of three.
%
% (2) 'GUI' → A logical true or false indicating whether or not the
% ScalarFit call originated from the GUI interface.
%
% (3) 'FittingParameters' → A vector of ones and zeros with the same
% size as the guess vector. In this vector, if the index value is equal to
% 1, then that index in the guess vector will be considered a free variable
% and will be fitted. If the index value is equal to 0, then that index in
% the guess vector will be frozen and not fitted.
% 
% (4) 'MaxFunEvals' → Specifies the maximum number of function evaluations
% in the call to fminsearch.
% 
% (5) 'OverridePlot' → Prevents MATLAB from plotting your data and fit.
% Should be entered in as a logical 1 or 0, or as "true" or "false".
% 
% (6) 'SizeNum' → Only required for the 'Rational' functional form.
% Specifies the number of terms in the numerator of the rational function.
%
% (7) 'SpecFcn' → Allows you to specify the function that you want to fit
% explicitly if it is not hard-coded into the script. You must choose
% 'Custom' as your functional form in order for this setting to work. The
% function specified must be a function of two vector variables c and x.
% Suppose you wanted to fit the function y = c1 + c2*x^2. The proper syntax
% for the function would be @(c,x) c(1) + c(2)*x.^2. It is generally
% appropriate to format the operations in the function handle as compatible
% with elementwise matrix operations.
% 
% FUNCTIONAL FORMS---------------------------------------------------------
% 
% Form - Guess - Function
% 
% 'Cubic' - [a,b,c,d] - ax^3+bx^2+cx+d
% 'Exponential' - [a,b,k] - ae^(bx)+k
% 'Gaussian' - [a,b,h,k] - ae^(-b(x-h)^2)+k
% 'Linear' - [a,b] - ax + b
% 'Logarithmic' - [a,h,k] - aln(x-h)+k
% 'Logistic' - [k,n,r] - nke^(rx)/(k-n+ne^(rx))
% 'Polynomial' - [a0,...,an] - sum(aix^i)
% 'Quadratic' - [a,h,k] - a(x-h)^2+k
% 'Rational' - [a0,...,am,b0,...,bn] - sum_0^m(aix^i)/sum_0^n(bix^i)
% 'Sinusoid' - [a,w,p,k] - asin(wx-p)+k

% /////////////////////////
% Benjamin Van Schaick
% Version 1.0.0
% Date Modified: 11/15/2024
% /////////////////////////

% Reformat to add c
% Add only fitting a subset of data in your data matrix.

% -----------------------(1) Parse Inputs----------------------------------

% (1.1) Data Matrix Manipulation

s = size(data);
if s(1) == 2 && s(2) ~= 2 % If true, data matrix is 2-by-n (n≠2)
    data = data'; % Transforms data matrix
    warning(sprintf(['Your data matrix is sideways (',num2str(s(1)),'-by-',num2str(s(2)),').\nYour data matrix has been transformed.']))
end

% (1.2) Ensure Guess Vector Integrity

if min(size(guess)) > 1 % If true, "guess" is a matrix, not a vector
    error(sprintf(['Your guess is ',num2str(s(1)),'-by-',num2str(s(2)),', and must be a vector.\nPlease input a vector corresponding to the chosen functional form.']))
end

if width(guess) == 1 && height(guess) ~= 1 % Guess vector is a column vector
    guess = guess'; % Transforms guess vector
    warning('Your guess vector was entered in as a column vector.\nIt has been transformed to a row vector')
end

% (1.3) Ensure Compatibility Between Guess and Form

l = length(guess);

if isequal(form,'Linear') || isequal(form,'LINEAR') || isequal(form,'linear')
    form = 'Linear'; % To collapse the three forms to one
    if width(guess) ~= 2 % "guess" doesn't have the required two elements for 'Linear'
        error(sprintf(['For the ',form,' functional form, your guess must have 2 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Quadratic') || isequal(form,'QUADRATIC') || isequal(form,'quadratic')
    form = 'Quadratic';
    if width(guess) ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Cubic') || isequal(form,'CUBIC') || isequal(form,'cubic')
    form = 'Cubic';
    if width(guess) ~= 4
        error(sprintf(['For the ',form,' functional form, your guess must have 4 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Polynomial') || isequal(form,'POLYNOMIAL') || isequal(form,'polynomial')
    form = 'Polynomial';
elseif isequal(form,'Rational') || isequal(form,'RATIONAL') || isequal(form,'rational')
    form = 'Rational';
elseif isequal(form,'Exponential') || isequal(form,'EXPONENTIAL') || isequal(form,'exponential')
    form = 'Exponential';
    if width(guess) ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Logarithmic') || isequal(form,'LOGARITHMIC') || isequal(form,'logarithmic')
    form = 'Logarithmic';
    if width(guess) ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Logistic') || isequal(form,'LOGISTIC') || isequal(form,'logistic')
    form = 'Logistic';
    if width(guess) ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Gaussian') || isequal(form,'GAUSSIAN') || isequal(form,'gaussian')
    form = 'Gaussian';
    if width(guess) ~= 4
        error(sprintf(['For the ',form,' functional form, your guess must have 4 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Sinusoid') || isequal(form,'SINUSOID') || isequal(form,'sinusoid')
    form = 'Sinusoid';
    if width(guess) ~= 4
        error(sprintf(['For the ',form,' functional form, your guess must have 4 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Custom') || isequal(form,'CUSTOM') || isequal(form,'custom')
    form = 'Custom';
end

% (1.4) Set Default Values

MaxFunEvals = 100000; % Maximum number of evaluations for fminsearch
fitgui = false; % Determines whether or not the ScalarFit call is from the GUI. If true, then ScalarFit will append the parameters argument with the error in the last index.
fitting_parameters = ones(size(guess));
override_plot = false;
specfcn = false;
bounds = cell(1,length(guess));
for i = 1:length(bounds)
    bounds{i} = [-Inf,Inf];
end

% (1.5) Assign Variable Values Based on Varargin

if nargin > 3
    if mod(nargin-3,2) == 0 || nargin == 4
        if nargin == 4
            varargin = varargin{1};
        end
        for i = 4:2:nargin
            if isequal(varargin{i-3},'Bounds')
                if isequal(class(varargin{i-2}),'char')
                    temp = str2num(varargin{i-2});
                else
                    temp = varargin{i-2};
                end
                for i2 = 1:3:length(temp)-2
                    bounds{temp(i2)} = [temp(i2+1) temp(i2+2)];
                end
            elseif isequal(varargin{i-3},'GUI')
                if isequal(varargin{i-2},'true') || isequal(varargin{i-2},'True') || isequal(varargin{i-2},'TRUE') || isequal(varargin{i-2},true)
                    fitgui = true;
                end
            elseif isequal(varargin{i-3},'FittingParameters')
                if isequal(class(varargin{i-2}),'char')
                    fitting_parameters = str2num(varargin{i-2});
                else
                    fitting_parameters = varargin{i-2};
                end
                if isequal(size(fitting_parameters),size(guess)) == false
                    error(sprintf('The size of your guess and fitting parameter vectors must be the same.\nMake sure that the fitting parameter vector is a row vector.'))
                end
            elseif isequal(varargin{i-3},'MaxFunEvals')
                if isequal(class(varargin{i-2}),'char')
                    MaxFunEvals = str2num(varargin{i-2});
                else
                    MaxFunEvals = varargin{i-2};
                end
            elseif isequal(varargin{i-3},'OverridePlot')
                if isequal(varargin{i-2},'true') || isequal(varargin{i-2},'True') || isequal(varargin{i-2},'TRUE') || isequal(varargin{i-2},true)
                    override_plot = true;
                end
            elseif isequal(varargin{i-3},'SizeNum')
                if isequal(class(varargin{i-2}),'char')
                    sizenum = str2num(varargin{i-2});
                else
                    sizenum = varargin{i-2};
                end
            elseif isequal(varargin{i-3},'SpecFcn')
                specfcn = varargin{i-2};
                if ~isequal(class(specfcn),'function_handle')
                    error(sprintf(['Your function for the SpecFcn name-value pair is a ',class(specfcn),' not a\n' ...
                        'function handle.']))
                end
                if ~isequal(form,'Custom')
                    error(sprintf(['You must specify the Custom functional form in order to use the SpecFcn\n' ...
                        'name-value pair. You chose the ',form,' functional form.']))
                end
            else
                error('Unsupported name in a name-value pair')
            end
        end
    end
end

options = optimset('MaxFunEvals',MaxFunEvals);

% ---------------------(2) Assign the Functional Form----------------------

if strcmp(form,'Cubic')
    fun = @(x) x(1)*(data(:,1)).^3 + x(2)*(data(:,1)).^2 + x(3)*(data(:,1)) + x(4);
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Custom')
    if isequal(specfcn,false)
        % MANUALLY HARD-CODE CUSTOM FUNCTION HERE IF YOU DON'T WANNA USE
        % 'SpecFcn' NAME-VALUE PAIR
        fun = @(x) x(1)*data(:,1)+x(2)+((x(3)*data(:,1)+x(4))./(x(5)*((data(:,1)).^2)+x(6)*data(:,1)+x(7)));
    else
        fun = @(x) feval(specfcn,x,data(:,1));
    end
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Exponential')
    fun = @(x) x(1)*exp(x(2)*data(:,1))+x(3);
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Gaussian')
    fun = @(x) x(1)*exp(-1*x(2)*(data(:,1)-x(3)).^2)+x(4);
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);
    disp('Note, the FWHM can be found using FWHM = -4*ln(2)/b')

elseif strcmp(form,'Linear')
    fun = @(x) x(1) * data(:,1) + x(2);
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Logarithmic')
    fun = @(x) x(1)*log(data(:,1)-x(2))+x(3);
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Logistic')
    fun = @(x) (x(2)*x(1)*exp(x(3)*data(:,1))).*((x(1)-x(2)+x(2)*exp(x(3)*data(:,1))).^(-1));
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Polynomial')
    fit_fun = @(x) polynomial(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Quadratic')
    fun = @(x) x(1) * (data(:,1) - x(2)).^2 + x(3);
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Rational')
    fit_fun = @(x) rational(guess.*(1-fitting_parameters)+x.*fitting_parameters);

elseif strcmp(form,'Sinusoid')
    if fitting_parameters(2) == 1 % If true, we are allowed to change the b value.
        if abs(max(diff(data(:,1))) - min(diff(data(:,1)))) < 1e-10 % Check to see if the independent variable is linearly spaced
            t_new = data(:,1);
            x_new = data(:,2);
        else
            t_new = linspace(min(data(:,1)),max(data(:,1)),1001)'; % new time vector with 1001 evenly spaced values
            x_new = interp1(data(:,1),data(:,2),t_new,'spline'); % interpolate x onto new time vector using spline method
        end
        dt = median(diff(t_new)); % time step
        Fs = 1/dt; % sampling frequency
        p = 2; % number of complex sinusoids in data
        f_est = rootmusic(x_new,p,Fs); % estimated frequencies

        % Uncomment to take a look at the pseudospectrum for visual verification of the frequency
%             [Pxx,f] = pmusic(x_new,p,[],Fs); % compute MUSIC pseudospectrum
%             plot(f,Pxx); % plot pseudospectrum
%             xlabel('Frequency (Hz)');
%             ylabel('Pseudospectrum');

        guess(2) = 2*pi*f_est(1);
    end

    fun = @(x) x(1)*sin(x(2)*data(:,1)-x(3))+x(4);
    fit_fun = @(x) fun(guess.*(1-fitting_parameters)+x.*fitting_parameters);

else
    error(sprintf('You did not enter a supported functional form. \nType "help ScalarFit" for more information.'))
end

p = fminsearch(@error_fit,guess,options);
p = guess.*(1-fitting_parameters)+p.*fitting_parameters;

% ------------------------(3) Function Dependencies------------------------

    function e = error_fit(p)
        e = sum((fit_fun(p) - data(:,2)).^2);
        for i3 = 1:length(bounds)
            temp = bounds{i3};
            if p(i3) < temp(1) || p(i3) > temp(2)
                e = Inf;
            end
        end
    end

    function y = polynomial(x)
        y = zeros(length(data(:,1)),1);
        n = length(x);
        for i1 = 0:n-1
            y = y + x(i1+1)*(data(:,1)).^i1;
        end
    end

    function y = rational(x)
        guess_num = x(1:sizenum);
        guess_denom = x(sizenum+1:length(guess));
        y = polynomial(guess_num)./polynomial(guess_denom);
    end

% ---------------------------(4) Plotting----------------------------------

if override_plot == false
    scatter(data(:,1),data(:,2),'b')
    hold on
    if isequal(form,'Rational')
        plot(linspace(min(data(:,1)),max(data(:,1)),10000),function_generator(linspace(min(data(:,1)),max(data(:,1)),10000),form,p,'SizeNum',sizenum),"Color",'r')
    elseif ~isequal(specfcn,false)
        plot(linspace(min(data(:,1)),max(data(:,1)),10000),function_generator(linspace(min(data(:,1)),max(data(:,1)),10000),form,p,'SpecFcn',specfcn),"Color",'r')
    else
        plot(linspace(min(data(:,1)),max(data(:,1)),10000),function_generator(linspace(min(data(:,1)),max(data(:,1)),10000),form,p),"Color",'r')
    end
    legend({'Data','Fit'},'location','best')
end

% ------------------------(5) Display Results------------------------------

disp(' ')
e = error_fit(p);
disp(['Error = ',num2str(e),' (SSE)'])
TSS = sum((mean(data(:,2)) - data(:,2)).^2);
r_squared = 1-(e/TSS);
disp(['R² = ',num2str(r_squared)])

if isequal(fitgui,'true')
    p = cat(2,p,e);
end

end