%%Computes the integration matrix including the sine function
%by Lisa Fischer
function [entry] = lagrange_sin_integrate( t, upper_bound, deltaT, c, k)

% degree of the Lagrange polynomial
n = length(t)-1;
entry = 0;

for l=n+1:-1:1    
    %compute sin(k*pi/deltaT *t(upper_bound)+r*pi/2)-sin(k*pi/deltaT *t(upper_bound-1)+r*pi/2)
    summandsSinPoly = zeros(1,l);
    
    for r=0:l-1
        %const independent of t
        const = 1/(k*pi/deltaT)^(r+1) * (factorial(l-1)/factorial(l-1-r));
        %function evaluation for t_{i}
        int_upper = t(upper_bound)^(l-1-r) *cos(k*pi/deltaT *t(upper_bound)+r*pi/2);
        %and t_{i-1}
        int_lower = t(upper_bound-1)^(l-1-r) *cos(k*pi/deltaT *t(upper_bound-1)+r*pi/2);
        %partial integrate summands:
        summandsSinPoly(r+1) = const *(int_upper  -  int_lower);       
    end
    
    %partial integration of polynomial and trigonometric function
    entry = entry - c(n+1-l+1)*sum(summandsSinPoly(:));
    
end