function frac = repetend2frac(str)
% Converts a decimal with repeating part into a fraction
% Input format examples:
% '0.12(34)'  -> 0.12343434...
% '12.(34)'   -> 12.343434...
% '0.(142857)'-> 0.142857142857...
%
% Output:
% frac = [numerator, denominator]

    % Find parentheses (repeating part)
    tokens = regexp(str,'^(\d*)\.?(\d*)\(?(\d*)\)?$','tokens');
    if isempty(tokens)
        error('Input format not recognized. Use form like "12.34(56)"');
    end
    tokens = tokens{1};
    intPart = tokens{1};       % integer part before decimal
    nonRep  = tokens{2};       % non-repeating decimals
    rep     = tokens{3};       % repeating decimals
    
    % Convert strings to numbers
    intPartNum = str2double(intPart);
    nonRepLen  = length(nonRep);
    repLen     = length(rep);
    
    % Case: no repeating part
    if repLen == 0
        numerator = intPartNum*10^nonRepLen + str2double(nonRep);
        denominator = 10^nonRepLen;
    else
        % Numerator = value without decimal up to repeating part minus value before repeating
        fullNum = str2double([nonRep rep]);
        if isempty(nonRep)
            preNum = 0;
        else
            preNum = str2double(nonRep);
        end
        numerator = fullNum - preNum + intPartNum*(10^(nonRepLen + repLen) - 10^nonRepLen);
        denominator = (10^repLen - 1) * 10^nonRepLen;
    end
    
    % Simplify fraction
    g = gcd(numerator, denominator);
    frac = [numerator/g, denominator/g];
end
