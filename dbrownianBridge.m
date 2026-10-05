%% derivative of Karhunen-Loève expansion
%by Lisa Fischer
%the first expansion term is not included here!
function [dbm] = dbrownianBridge(deltaT, t, xi)

dbm =zeros(size(xi,1),1);
for k=1:size(xi,2)
    dbm = dbm + cos(k*pi*t/deltaT)*xi(:,k);%fix this for multidim sde
end
%dbm = sqrt(2/deltaT) .* dbm;

end