%% coeff
%this function calculates the coefficients to integrate the lagrange
%polynomial
%out of these we will later on get one entry of the integration matrix
%by Felix Binkowski


function [ coefficients ] = coeff( t, nodes, ref_index)
% ref_index is the index of the Langrange polynomial

% save t for integration
t_nodes=t;

% nodes represented by the quadrature nodes used
if strcmp(nodes, 'legendre')
    t_nodes(1)=[];
    t_nodes(end)=[];
elseif strcmp(nodes, 'radau')
    t_nodes(1)=[];
end


%% compute constants for lagrange polynomial
c=zeros(1,length(t_nodes));
c(end-1)=1;
c(end)=1;

for i=1:length(t_nodes)
    if i~=ref_index
        c(end-1)=c(end-1)*(t_nodes(ref_index)-t_nodes(i));
        c(end)=c(end)*t_nodes(i);
    end
end

if(length(t_nodes) >2)
    t_tmp=t_nodes;
    t_tmp(ref_index)=[];
    
    for j=0:length(t_tmp)-1
        c(1) = c(1) + t_tmp(j+1);
        if j==0
            for i=2:length(t_tmp)-1
                c(i)= c(i)+c(i-1)*t_tmp(i+j);
            end
        else
            for i=2:length(t_tmp)-j
                c(i)= c(i)+c(i-1)*t_tmp(i+j);
            end
        end
    end
end
%create coeeficient vector with all entries

c_final=zeros(1,length(t_nodes));
c_final(1)=c(end);
if(length(t_nodes)>2)
    c_final(2:end-1)=fliplr(c(1:end-2));
end
c_final(end)=1;
c_final(end+1)=c(end-1);

coefficients = c_final;

end