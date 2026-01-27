close all
clear all

dt = 1/4; % length of time step
dW =  randn * sqrt(dt); % Wiener process increment

nBridges = 3;  % number of different bridges to plot
M = 4; % number of terms in the bridge (M=1 only linear interpolation)

bridgeRange=1:M;

xi = randn(nBridges, M-1);

subdt = 1/128;

t=0:subdt:dt;

% [nodes, ~] = lgwt(5,0,dt);
% tc = nodes(4);
% 
% v(1) = dW*tc/sqrt(dt);
% for b=2:M
%     v(b) = brownianBridge(dW, dt, tc, xi(1,1:b-1));
% end

for b=bridgeRange

figure(b)
hold on
plot(t, dW*t/sqrt(dt)); % linear only

for i=1:nBridges
    bridge = brownianBridge(dW, dt, t, xi(i,1:b-1));
    plot(t, bridge);
end

end

shift = @(t) 2*t/dt-1;
J{1} = @(t) 0.25*(shift(t).^2-1);
J{2} = @(t) 0.5*(shift(t).^3-shift(t));
J{3} = @(t) 15*shift(t).^4/16 - 9*shift(t).^2/8 + 3/16;
J{4} = @(t) 7*shift(t).^5/4 - 5*shift(t).^3/2 + 3*shift(t)/4;

for b=1:4
    figure(b+4)
    hold on
    plot(t, dW*t/sqrt(dt))
    
    for i=1:nBridges
        bridge = dW*t/sqrt(dt);
        for bb=1:b-1
            bridge = bridge + sqrt(1/((bb+1)*(bb+2)))*xi(i,bb)*J{bb}(t);
        end
        plot(t, bridge);
    end
end