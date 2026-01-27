%% lagrange_integrate
%calculates one spectral integration matrix entry
%via integration over the lagrange polynomial and according coefficients
%by Felix Binkowski

function [entries] = lagrange_integrate( t, upper_bound, c, nodes)

%matrix size depends on nodes
if strcmp(nodes, 'radau')
    len = length(t)-1;    
elseif strcmp(nodes, 'legendre')
    len = length(t)-2;    
else
    len = length(t);
end

entries = 0;
    
%calculate lagrange integration     
switch mod(len, 2)
    case 0
        
        for i=1:len
            entries = entries + ((-1)^(i)*(1/i)) *c(i) *(t(upper_bound)^(i)-t(upper_bound-1)^(i));
        end
        
        
    otherwise
        for i=1:len
            entries = entries + ((-1)^(i+1)*(1/i)) *c(i) *(t(upper_bound)^(i)-t(upper_bound-1)^(i));
        end
        
        
end
          
entries = entries/c(end);

end
