%% derivative of Karhunen-Loève expansion
%by Lisa Fischer
%the first expansion term is not included here!
function [dbm] = dbrownianBridge(deltaT, t, xi)

sumBm =zeros(size(xi,1),1);
for k=1:length(xi)
    sumBm = sumBm + cos(k*pi*t/deltaT)*xi(:,k);%fix this for multidim sde
end
dbm = sqrt(2/deltaT) .*sumBm;

end