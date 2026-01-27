function [S, t, x] = computeSpecMat(t_begin, t_end, step_size, intervals, col_points, nodes)

% compute spectral integration matrix
if (col_points ==1)
    disp('the method needs at least 2 collocation points');
else
    t = col_nodes(t_begin, t_begin+step_size, col_points, nodes);
    S = spectral_matrix(t, nodes);
end

if (strcmp (nodes, 'radau'))
    x = linspace(t_begin, t_end, (col_points)*intervals+1);
elseif strcmp (nodes, 'legendre')
    x = linspace(t_begin, t_end, (col_points+1)*intervals+1);
else
    x = linspace(t_begin, t_end, (col_points-1)*intervals+1);
end

end