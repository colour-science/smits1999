function [f, g]  = colorFun(c)
 
global A
global x
global An

y = x + An * c;

if(max(y) > 1.0)  % penalty for values greater than 1
   mx = 100 * (max(y) - 1);
else 
   mx = 0;
end

f = 10 * norm(diff(y)) + mx;   % function to minimize

g = [-00-min(y)  max(y)-1.2];  % constrain the spectrum between 0 and 1.2
