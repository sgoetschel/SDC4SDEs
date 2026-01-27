%% spectral integration matrix
%computation of teh spectral integration matrix

function [ S ] = spectral_matrix( t, nodes)

len = length(t);

S = zeros(len-1, len);

%spectral integration matrix w.r.t time grid
if strcmp(nodes, 'radau')   
    
       for col=2:len
            c = coeff(t, nodes, col-1);
          for row=1:len-1 
            S(row,col) = lagrange_integrate( t, row+1, c, nodes);
          end
       end 
     
    
elseif strcmp(nodes, 'legendre')
    
       for col=2:len-1
           c = coeff(t, nodes, col-1);
         for row=1:len-1
            S(row,col) = lagrange_integrate( t, row+1, c, nodes);
         end
       end
       
else
       for col=1:len
           c = coeff(t, nodes, col);
        for row=1:len-1  
           S(row,col) = lagrange_integrate( t, row+1, c, nodes);
        end
       end 
end



end

