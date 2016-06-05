function res = cinput(s,DefaultValue)
% 
% res = cinput(s,DefaultValue)
% User input. If no value is given, res takes the DefaultValue.  

text = [s  ' (default '  num2str(DefaultValue)  ')= ']; 
res = input(text); 
if isempty(res)
    res = DefaultValue; 
end

