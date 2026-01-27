%% Karhunen-Loève expansion - Browninan bridge
% by Lisa Fischer
%the first expansion term is not included here!
function [bm] = brownianBridge(eta, deltaT , t, xi)

sum = 0;
for k=1:length(xi)
    sum = sum + sin(k*pi*t/deltaT)/k*xi(k);
end
bm = eta*t/sqrt(deltaT) + sqrt(2*deltaT)/pi .*sum;
%bm = eta*t + sqrt(2*deltaT)/pi .*sum; am 09.07. geändert

end