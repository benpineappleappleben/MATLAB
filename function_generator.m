function y = function_generator(domain,form,parameters,varargin)
% DESCRIPTION--------------------------------------------------------------
% 
% function_generator that generates a plot of a function given a functional
% form as well as the parameters necessary to uniquely define that
% function.
% 
% ARGUMENTS----------------------------------------------------------------
% 
% The domain argument  must be a row vector. If it is entered as a column
% vector, it will be transformed into a row vector. Consider using linspace
% in generating your domain for you. A linearly spaced domain is not
% required.
% 
% The first argument of FIT1D is the functional form, which must be entered
% as a character string. This argument chooses the function that will be
% used to generate a trace. A list of all functional forms can be found at
% the end of this help article as well as in the complementary document
% called "FIT1D Functional Forms v1.0.0.pdf".
% 
% The parameters argument must be a row or column vector containing all of
% the parameters necessary to generate the function of the desired form. A
% list of functional forms and their corresponding parameters can be found
% in the document entitled "FIT1D Functional Forms v1.0.0.pdf".
%
% VARARGIN ARGUMENTS-------------------------------------------------------
% 
% All varargin inputs should be entered in a name-value format as follows:
% function_generator(form,data,guess,'Name',value,'Name',value,...)
%
% Below is a list of the supported name-value pairs along with what they
% change in the function_generator pipeline. function_generator does not
% require a certain number of varargin arguments to operate unless you are
% using a form that requires a name-value pair. For example, the 'Rational'
% functional form requires a name-value pair to operate.
%
% In alphabetical order, the supported name-value pairs are:
%
% (1) 'SizeNum' → Only required for the 'Rational' functional form.
% Specifies the number of terms in the numerator of the rational function.
%
% (2) 'SpecFcn' → Allows you to specify the function that you want to fit
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

% -------------------------(1) Parse Inputs--------------------------------

% (1.1) Ensure Integrity of Domain Argument

s = size(domain);
if min(size(domain)) > 1
    error(['Your domain is ',num2str(s(1)),'-by-',num2str(s(2)),'. The domain must be a vector.'])
elseif height(domain) > width(domain)
    domain = domain';
    warning('Your domain must be a row vector, not a column vector. It has been transformed.')
end

% (1.2) Ensure Integrity of Parameters Argument

p = parameters;
if min(size(p)) > 1
    error(['Your parameters is ',num2str(s(1)),'-by-',num2str(s(2)),'. The domain must be a vector.'])
elseif height(p) > width(p)
    p = p';
    warning('Your parameters must be a row vector, not a column vector. It has been transformed.')
end

% (1.3) Ensure Compatibility Between Guess and Form

l = width(p);

if isequal(form,'Linear') || isequal(form,'LINEAR') || isequal(form,'linear')
    form = 'Linear'; % To collapse the three forms to one
    if l ~= 2 % "guess" doesn't have the required two elements for 'Linear'
        error(sprintf(['For the ',form,' functional form, your guess must have 2 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Quadratic') || isequal(form,'QUADRATIC') || isequal(form,'quadratic')
    form = 'Quadratic';
    if l ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Cubic') || isequal(form,'CUBIC') || isequal(form,'cubic')
    form = 'Cubic';
    if l ~= 4
        error(sprintf(['For the ',form,' functional form, your guess must have 4 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Polynomial') || isequal(form,'POLYNOMIAL') || isequal(form,'polynomial')
    form = 'Polynomial';
elseif isequal(form,'Rational') || isequal(form,'RATIONAL') || isequal(form,'rational')
    form = 'Rational';
elseif isequal(form,'Exponential') || isequal(form,'EXPONENTIAL') || isequal(form,'exponential')
    form = 'Exponential';
    if l ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Logarithmic') || isequal(form,'LOGARITHMIC') || isequal(form,'logarithmic')
    form = 'Logarithmic';
    if l ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Logistic') || isequal(form,'LOGISTIC') || isequal(form,'logistic')
    form = 'Logistic';
    if l ~= 3
        error(sprintf(['For the ',form,' functional form, your guess must have 3 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Gaussian') || isequal(form,'GAUSSIAN') || isequal(form,'gaussian')
    form = 'Gaussian';
    if l ~= 4
        error(sprintf(['For the ',form,' functional form, your guess must have 4 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Sinusoid') || isequal(form,'SINUSOID') || isequal(form,'sinusoid')
    form = 'Sinusoid';
    if l ~= 4
        error(sprintf(['For the ',form,' functional form, your guess must have 4 elements.\nYour guess has ',num2str(l),' elements.']))
    end
elseif isequal(form,'Custom') || isequal(form,'CUSTOM') || isequal(form,'custom')
    form = 'Custom';
end

% (1.4) Assign Variable Values Based on Varargin

specfcn = false;

if nargin > 3
    if mod(nargin-3,2) == 0
        for i = 4:2:nargin
            if isequal(varargin{i-3},'SizeNum')
                sizenum = varargin{i-2};
                sizeden = length(parameters) - sizenum;
            elseif isequal(varargin{i-3},'SpecFcn')
                specfcn = varargin{i-2};
            else
                error('Unsupported name in a name-value pair')
            end
        end
    else
        error('Unexpected number of inputs. Each name and value should be inputted as a pair')
    end
end

if strcmp(form,'Linear')
    y = p(1) * domain + p(2) * ones(s);

elseif strcmp(form,'Quadratic')
    y = ((domain - p(2) * ones(s)).^2) * p(1) + p(3) * ones(s);

elseif strcmp(form,'Cubic')
    y = p(1) * domain.^3 + p(2) * domain.^2 + p(3) * domain + p(4) * ones(s);

elseif strcmp(form,'Polynomial')
    n = length(p);
    y = zeros(1,length(domain));
    for i = 0:n-1
        y = y + p(i+1)*(domain).^i;
    end
elseif strcmp(form,'Rational')
    y_num = zeros(1,length(domain));
    y_denom = y_num;
    p_denom = p(1+sizenum:length(p));
    for i = 0:sizenum-1
        y_num = y_num + p(i+1)*(domain).^i;
    end
    for i = 0:sizeden-1
        y_denom = y_denom + p_denom(i+1)*(domain).^i;
    end
    y = y_num./y_denom;
elseif strcmp(form,'Exponential')
    y = p(1) * exp(p(2) * domain) + p(3) * ones(s);

elseif strcmp(form,'Logarithmic')
    y = p(1)*log(domain-p(2))+p(3);

elseif strcmp(form,'Logistic')
    y = (p(2)*p(1)*exp(p(3)*domain)).*((p(1)-p(2)+p(2)*exp(p(3)*domain)).^(-1));

elseif strcmp(form,'Gaussian')
    y = p(1) * exp(-1 * p(2) * (domain - p(3) * ones(s)).^2) + p(4) * ones(s);
    
elseif strcmp(form,'Sinusoid')
    y = p(1) * sin(p(2) * domain - p(3)) + p(4);

elseif strcmp(form,'Custom')
    if isequal(specfcn,false)
        % HARD-CODE CUSTOM FUNCTION HERE AS WELL FOR PLOTTING (see FIT1D)
        y = p(1)*domain+p(2)+((p(3)*domain+p(4))./(p(5)*((domain).^2)+p(6)*domain+p(7)));
    else
        y = feval(specfcn,p,domain);
    end
else
    error(sprintf(['You did not enter a supported functional form.\nType "help function_generator" for more information']))
end
end
