%% col nodes
% create collocation nodes
% by Felix Binkowski

function [t] = col_nodes(a,b,t_steps,nodes)

if strcmp(nodes,'legendre')
    % Gauss-Legendre nodes
    [t(1:t_steps),~]=lgwt(t_steps,a,b);
    t_tmp(1:t_steps)=flipud(t(1:t_steps)');
    %t_tmp = [a,t(1:t_steps),b]
    t(1:t_steps+2)=[a,t_tmp,b];
    
elseif strcmp(nodes,'lobatto')
    % Lobatto nodes
    % A-stable (left half plane) and not L-stable => not good for
    % stiff problems
    [t,~,~]=lglnodes(t_steps-1,a,b);
    t=flipud(t)';
    
elseif strcmp(nodes,'radau')
    % Radau nodes
    % L-stable => stiff problems
    % ( A-stable and lim |Re(z)| = 0 with Re(z) -> inf )
    % includes right interval border
    [t(1:t_steps),~,~]=lgrnodes(t_steps-1,a,b);
    t(t_steps+1) = b;
    
elseif strcmp(nodes,'radau_left')
    % Radau nodes
    % L-stable => stiff problems
    % ( A-stable and lim |Re(z)| = 0 with Re(z) -> inf )
    % includes left interval border
    [t(1:t_steps),~,~]=lgrnodes(t_steps-1,a,b);
    %t(t_steps+1) = b;
    %t_tmp = flipud(t(1:t_steps));
    t = [a , t(2:t_steps)];
    
elseif strcmp(nodes,'equidist')
    % equidistant time grid
    t=linspace(a,b,t_steps);
end

end

