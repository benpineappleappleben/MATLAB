function y = unit2unit(value,init_unit,fin_unit)
% unit2unit takes an initial value with units init_unit (string) and
% performs the necessary dimensional analysis to convert it into a new
% value with units fin_unit (string). Value is allowed to be a vector in
% the case that you want to convert all elements of that vector from one
% unit to another.

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Created by: Benjamin Van Schaick
% Date Modified: 4/22/2026
% Version: 1.0.6
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Physical constants defined by SI
h = 6.62607015*10^(-34); % Planck's constant, exactly, in J*s
c = 299792458; % Speed of light, exactly, in m/s

% Parse unit strings for simplicity later
init_unit = unitParse(init_unit);
fin_unit = unitParse(fin_unit);

% Perform the desired transformation
if isequal(init_unit,'eV') && isequal(fin_unit,'J')
    y = 1.60218e-19*value;
elseif isequal(init_unit,'eV') && isequal(fin_unit,'nm')
    y = h*c*10^9./(value*1.60218e-19);
elseif isequal(init_unit,'eV') && isequal(fin_unit,'wavenumber')
    y = value*1.60218e-19/(h*100*c);
elseif isequal(init_unit,'Gauss') && isequal(fin_unit,'T')
    y = value/10000;
elseif isequal(init_unit,'Hertz') && isequal(fin_unit,'wavenumber')
    y = (1/(100*c))*value;
elseif isequal(init_unit,'Hex') && isequal(fin_unit,'RGB')
    if isequal(class(value),'char')
        if isequal(value(1),'#')
            value(1) = [];
        end
    elseif isequal(class(value),'string')
        value = char(value);
        if isequal(value(1),'#')
            value(1) = [];
        end
    elseif isequal(class(value),'double')
        value = num2str(value);
    else
        error(['The class of your hexadecimal color code must be character, string, or double. The class of your hex code is ',class(value),'.'])
    end
    value = [hex2dec(value(1:2)),hex2dec(value(3:4)),hex2dec(value(5:6))];
    y = (1/255)*value;
elseif isequal(init_unit,'J') && isequal(fin_unit,'eV')
    y = value/1.60218e-19;
elseif isequal(init_unit,'J') && isequal(fin_unit,'wavenumber')
    y = (1/(100*h*c))*value;
elseif isequal(init_unit,'nm') && isequal(fin_unit,'eV')
    y = h*c*10^9./(value*1.60218e-19);
elseif isequal(init_unit,'nm') && isequal(fin_unit,'wavenumber')
    y = 10^7./value;
elseif isequal(init_unit,'rad/s') && isequal(fin_unit,'wavenumber')
    y = (1/(200*pi*c))*value;
elseif isequal(init_unit,'RGB') && isequal(fin_unit,'Hex')
    value = round(255*value);
    y = [dec2hex(value(1),2),dec2hex(value(2),2),dec2hex(value(3),2)];
elseif isequal(init_unit,'T') && isequal(fin_unit,'Gauss')
    y = value*10000;
elseif isequal(init_unit,'wavenumber') && isequal(fin_unit,'eV')
    y = h*100*c*value/1.60218e-19;
elseif isequal(init_unit,'wavenumber') && isequal(fin_unit,'Hertz')
    y = 100*c*value;
elseif isequal(init_unit,'wavenumber') && isequal(fin_unit,'J')
    y = 100*h*c*value;
elseif isequal(init_unit,'wavenumber') && isequal(fin_unit,'nm')
    y = 10^7./value;
elseif isequal(init_unit,'wavenumber') && isequal(fin_unit,'nm^(-1)')
    y = value*10e-7;
elseif isequal(init_unit,'wavenumber') && isequal(fin_unit,'rad/s')
    y = 200*pi*c*value;
else
    error(['Unit conversion from ',init_unit,' to ',fin_unit,' is not supported.'])
end

function v = unitParse(u)
    token_warning = 0;
    v = u;
    if isequal(u,'eV') || isequal(u,'electronvolts') || isequal(u,'electronvolt') || isequal(u,'ElectronVolts') || isequal(u,'ElectronVolt') || isequal(u,'ev')
        if ~isequal(u,'eV')
            v = 'eV';
            token_warning = 'energy in eV';
            token_standard = v;
        end
    elseif isequal(u,'Gauss') || isequal(u,'gauss')
        if ~isequal(u,'Gauss')
            v = 'Gauss';
            token_warning = 'magnetic field strength in Gauss';
            token_standard = v;
        end
    elseif isequal(u,'Hertz') || isequal(u,'hertz') || isequal(u,'s^-1') || isequal(u,'s^(-1)') || isequal(u,'s^{-1}')
        if ~isequal(u,'Hertz')
            v = 'Hertz';
            token_warning = 'frequency in Hertz';
            token_standard = v;
        end
    elseif isequal(u,'Hex') || isequal(u,'hex') || isequal(u,'Hex Code') || isequal(u,'HexCode') || isequal(u,'Hexadecimal') || isequal(u,'hexadecimal')
        if ~isequal(u,'Hex')
            v = 'Hex';
            token_warning = 'color code in hexadecimal';
            token_standard = v;
        end
    elseif isequal(u,'J') || isequal(u,'Joules') || isequal(u,'joules')
        if ~isequal(u,'J')
            v = 'J';
            token_warning = 'energy in Joules';
            token_standard = v;
        end
    elseif isequal(u,'micron') || isequal(u,'microns') || isequal(u,'um')
        if ~isequal(u,'um')
            v = 'um';
            token_warning = 'wavelength in micrometers';
            token_standard = v;
        end
    elseif isequal(u,'nm^(-1)') || isequal(u,'nm^-1')
        if ~isequal(u,'nm^(-1)')
            v = 'nm^(-1)';
            token_warning = 'wavenumber in inverse nanometers';
            token_standard = v;
        end
    elseif isequal(u,'omega') || isequal(u,'rad/s') || isequal(u,'rad/sec') || isequal(u,'angular frequency')
        if ~isequal(u,'rad/s')
            v = 'rad/s';
            token_warning = 'angular frequency in radians per second';
            token_standard = v;
        end
    elseif isequal(u,'RGB') || isequal(u,'rgb')
        if ~isequal(u,'RGB')
            v = 'RGB';
            token_warning = 'color code in RGB';
            token_standard = v;
        end
    elseif isequal(u,'Tesla') || isequal(u,'tesla') || isequal(u,'T')
        if ~isequal(u,'T')
            v = 'T';
            token_warning = 'magnetic field strength in Tesla';
            token_standard = v;
        end
    elseif isequal(u,'Wavelength') || isequal(u,'wavelength') || isequal(u,'nm')
        if ~isequal(u,'nm')
            v = 'nm';
            token_warning = 'wavelength in nanometers';
            token_standard = v;
        end
    elseif isequal(u,'Wavenumber') || isequal(u,'wavenumber') || isequal(u,'cm^-1') || isequal(u,'cm^(-1)') || isequal(u,'cm^{-1}')
        if ~isequal(u,'wavenumber')
            v = 'wavenumber';
            token_warning = 'wavenumber';
            token_standard = v;
        end
    else
        error(['Unit "',u,'" not recognized. Verify that your units are supported.'])
    end
    if ~isequal(token_warning,0)
        warning(['Your unit (',u,') was taken to mean ',token_warning,' and was converted to the system standard (',token_standard,').'])
    end
end
end