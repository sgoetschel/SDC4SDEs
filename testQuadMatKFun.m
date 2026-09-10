%% testQuadMatKFun.m
% Correctness tests for Q: a matrix whose (i,j) entry is claimed to be
%
%       Q(i,j) = integral_{t_i}^{t_{i+1}}  ell_j(t) * cos(k*pi*t/L) dt
%
% where ell_j is the j-th Lagrange basis polynomial on the GLOBAL node
% set t_1..t_N, and L = t1 - t0 is the length of the whole domain.
%
% STRATEGY
% -----------------------------------------------------------------------
% Rather than only checking exactness on a handful of basis functions
% (as in test_computeQkcMat.m), this script reconstructs EVERY entry of
% Q independently:
%   1. Evaluate ell_j via the barycentric Lagrange formula (numerically
%      stable for any node count/spacing -- avoids the ill-conditioning
%      of an explicit monomial expansion).
%   2. Integrate ell_j(t)*cos(k*pi*t/L) over each segment [t_i,t_{i+1}]
%      using MATLAB's adaptive `integral` (the integrand is smooth/
%      entire, so this converges to near machine precision).
%   3. Compare entrywise against the actual Q under test.
% This is a ground-truth structural check, independent of however Q was
% actually constructed (monomial expansion, recursion, series, etc).
%
% ASSUMPTIONS -- EDIT AS NEEDED
% -----------------------------------------------------------------------
% - Calling convention: Q = computeQMat(t0, t1, step_size, intervals, nodes, k),
%   with node locations optionally returned as a 2nd output. Edit the
%   CALL FUNCTION block if your actual function name/signature differs.
% - NODE_KIND used only as a fallback if node locations aren't returned.

clear; clc; close all;

%% ------------------------- CONFIG --------------------------------------
t0        = 0;
t1        = 0.1; %2*pi;
intervals = 1;
step_size = (t1 - t0)/intervals;
nnodes     = 4;                
nodeType = 'lobatto';           

L           = t1 - t0;
K_TEST_LIST = [1, 2, 3, 4, 5, 10, 50];    % <-- EDIT: which k values to check
Kdesign = K_TEST_LIST(end);

tol_entry  = 1e-9;    % tolerance for entrywise reconstruction
tol_smooth = 1e-6;    % tolerance for the end-to-end functional check

%% ------------------------- BUILD MATRIX --------------------------------
t = col_nodes(t0, t0+step_size, nnodes, nodeType);
Qkc = quadMatKFun(t, step_size, nodeType, Kdesign);

N = size(Qkc,2);  % for each k (i.e., frequency), Qkc is (N-1)xN, where N = nnodes

t     = t(:);
tprev = t(1:end-1);
tnext = t(2:end);    
w     = barycentricWeights(t);

%% ==================== TEST 0: partition-of-unity row check ================ 
fprintf('\n============================== Test 0: row sums (partition of unity) ==============================\n');

for k = K_TEST_LIST
    fprintf('k = %d:\n', k);



    % sum_j ell_j(t) == 1 identically, so sum_j Q(i,j) MUST equal
    % integral_{t_i}^{t_{i+1}} cos(k*pi*t/L) dt. Cheap necessary condition,
    % checked before the expensive full entrywise reconstruction.
    rowsum_approx = sum(Qkc(:,:,k),2);
    rowsum_exact  = cosIntegral(k, tprev, tnext, L);
    err0 = max(abs(rowsum_approx - rowsum_exact));
    fprintf('  max abs row-sum error = %.3e  [%s]\n', err0, tern(err0<tol_entry,'PASS','FAIL'));
end

%% ==================== TEST 1: full entrywise reconstruction ================
fprintf('\n==============================  Test 1: entrywise reconstruction (barycentric ell_j + quadrature) ==============================\n');
    
for k = K_TEST_LIST
    fprintf('k = %d:\n', k);

    Qref = zeros(N-1,N);
    for j = 1:N
        integrand = @(x) evalLagrangeBasis(j, x, t, w) .* cos(k*pi*x/L);
        for i = 1:N-1
            Qref(i,j) = integral(integrand, tprev(i), tnext(i), 'AbsTol', 1e-13, 'RelTol', 1e-13);
        end
    end
    E = abs(Qkc(:,:,k) - Qref);
    [maxErr, linIdx] = max(E(:));
    [iWorst, jWorst] = ind2sub(size(E), linIdx);
    fprintf('  max abs entry error   = %.3e  at (i,j) = (%d,%d)  [%s]\n', ...
        maxErr, iWorst, jWorst, tern(maxErr<tol_entry,'PASS','FAIL'));
    fprintf('  mean abs entry error  = %.3e\n', mean(E(:)));
    
    % figure;
    % imagesc(log10(max(E,eps))); colorbar;
    % xlabel('j (basis index)'); ylabel('i (segment index)');
    % title(sprintf('log10(|Q - Q_{ref}|), k = %d', k));
end

%% ==================== TEST 2: end-to-end functional check ==================
fprintf('\n==============================  Test 2: end-to-end check on sampled functions ==============================\n');
for k = K_TEST_LIST
    fprintf('k = %d:\n', k);
    % Sanity check using actual sampled functions (not just isolated basis
    % columns): Q*g(nodes) should match a direct high-accuracy quadrature of
    % g(t)*cos(k*pi*t/L), for g NOT necessarily equal to the cosine itself.
    testFcns = {@(x) exp(0.2*x), @(x) sin(3*x+0.4), @(x) x.^2 - x + 1};
    names    = {'exp(0.2t)', 'sin(3t+0.4)', 't^2-t+1'};
    for m = 1:numel(testFcns)
        fh = testFcns{m};
        g  = fh(t);
        Iapprox = Qkc(:,:,k)*g;
        Iexact  = zeros(N-1,1);
        for i = 1:N-1
            Iexact(i) = integral(@(x) fh(x).*cos(k*pi*x/L), tprev(i), tnext(i), ...
                                  'AbsTol', 1e-12, 'RelTol', 1e-12);
        end
        err = max(abs(Iapprox - Iexact));
        fprintf('  %-14s max segment error = %.3e  [%s]\n', names{m}, err, tern(err<tol_smooth,'PASS','FAIL'));
    end
end % k loop

%% ------------------------- helper functions ---------------------------------
function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end

function I = cosIntegral(k, a, b, L)
    % integral_a^b cos(k*pi*t/L) dt, elementwise for vector a,b
    if k == 0
        I = b - a;
    else
        omega = k*pi/L;
        I = (sin(omega*b) - sin(omega*a)) / omega;
    end
end

function w = barycentricWeights(t)
    n = numel(t);
    w = ones(n,1);
    for j = 1:n
        w(j) = 1/prod(t(j) - t([1:j-1, j+1:n]));
    end
end

function y = evalLagrangeBasis(j, x, t, w)
    % Evaluate the j-th Lagrange basis polynomial at points x via the
    % (numerically stable) barycentric formula:
    %   ell_j(x) = (w_j/(x-t_j)) / sum_i (w_i/(x-t_i))
    % with exact handling at x == t_m (Kronecker delta), needed because
    % the formula has a removable singularity there.
    x = x(:)';
    y = zeros(size(x));
    tol = 10*eps(max(abs(t)));
    for idx = 1:numel(x)
        diffs = x(idx) - t;
        hit = find(abs(diffs) < tol, 1);
        if ~isempty(hit)
            y(idx) = (hit == j);
        else
            terms = w ./ diffs;
            y(idx) = terms(j) / sum(terms);
        end
    end
end

function t = getNodeLocations(t0, step_size, intervals, nodes, nodeKind)
    if isnumeric(nodes) && isscalar(nodes)
        switch nodeKind
            case 'LGL',        ref = lglnodes(nodes);
            case 'CGL',        ref = cglnodes(nodes);
            case 'equispaced', ref = linspace(-1,1,nodes)';
            otherwise, error('Unknown NODE_KIND: %s', nodeKind);
        end
    elseif isnumeric(nodes) && isvector(nodes)
        ref = nodes(:);
    else
        error('Edit getNodeLocations to handle nodes of class %s', class(nodes));
    end
    ref = sort(ref);
    t = [];
    for k = 0:intervals-1
        a  = t0 + k*step_size;
        b  = a + step_size;
        tk = a + (b-a)*(ref+1)/2;
        if k > 0
            tk(1) = [];
        end
        t = [t; tk]; %#ok<AGROW>
    end
end

function x = lglnodes(n)
    if n == 1, x = 0; return; end
    N1 = n-1;
    x  = cos(pi*(0:N1)/N1)';
    P  = zeros(n, n);
    xold = 2*ones(size(x));
    while max(abs(x-xold)) > eps
        xold = x;
        P(:,1) = 1; P(:,2) = x;
        for kk = 2:N1
            P(:,kk+1) = ((2*kk-1)*x.*P(:,kk) - (kk-1)*P(:,kk-1))/kk;
        end
        x = xold - (x.*P(:,n) - P(:,N1)) ./ (n*P(:,n));
    end
    x = sort(x);
end

function x = cglnodes(n)
    if n == 1, x = 0; return; end
    j = (0:n-1)';
    x = sort(cos(pi*j/(n-1)));
end