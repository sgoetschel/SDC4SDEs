%% testSpecMat.m
% Correctness tests for computeSpecMat(t0, t1, step_size, intervals, nodes)
%
% WHAT THIS TESTS
% -----------------------------------------------------------------------
% A spectral/collocation integration matrix S is assumed to map a vector
% of function samples F = f(t_1),...,f(t_N) taken at the collocation
% nodes to the running (cumulative) integral:
%
%       I(i) = integral from t0 to t_i of f(tau) dtau ,   i = 1..N
%       I ≈ S * F
%
% This is the standard convention for "integration matrices" in
% pseudospectral / spectral collocation methods (the counterpart of a
% differentiation matrix D, for which D*F ≈ f'). If your matrix instead
% returns, e.g., the integral over each sub-interval only (rather than
% cumulative from t0), adjust the "exact value" computations in Tests
% 1/3/4 accordingly -- the rest of the harness (calling the function,
% getting node locations, reporting errors) stays the same.
%
% HOW TO USE
% -----------------------------------------------------------------------
%   1. Edit the CONFIG block below to match a typical call you would
%      make to your own function.
%   2. If computeSpecMat does NOT return the node locations as a second
%      output, edit getNodeLocations() (bottom of this file) so that it
%      reconstructs t_1..t_N exactly the way your function places them.
%      The default assumes `nodes` is either an integer node count
%      (defaulting to Legendre-Gauss-Lobatto points) or an explicit
%      vector of reference nodes on [-1,1], affinely mapped onto each
%      sub-interval.
%   3. Run the script.
%
% All tests print PASS/FAIL against a tolerance plus the actual error,
% so you can judge borderline cases yourself rather than trusting a
% single boolean.

clear; clc; close all;

%% ------------------------- CONFIG --------------------------------------
t0        = 0;
t1        = 0.5; % 2*pi;
intervals = 1;                        % number of sub-intervals
step_size = (t1 - t0)/intervals;      % keep consistent with t0,t1,intervals
nnodes     = 4;                        % <-- EDIT: whatever `nodes` means for
                                       %     your function (a count, a
                                       %     string like 'LGL', or an
                                       %     explicit vector of reference
                                       %     points on [-1,1])
nodeType = 'lobatto';

SHARED_BOUNDARY_NODES = true;         % true if consecutive sub-intervals
                                       % share their endpoint node
                                       % (typical for Lobatto-type
                                       % schemes; set false for Gauss-type
                                       % schemes with no endpoints)

tol_poly   = 1e-9;   % tolerance for exact polynomial integration
tol_smooth = 1e-6;   % tolerance for smooth transcendental functions

%% ------------------------- CALL FUNCTION --------------------------------
[S, t] = computeSpecMat(t0, t1, step_size, intervals, nnodes, nodeType);

N = size(S,2); % S is (N-1)xN, N = nnodes

%% ------------------------- NODE LOCATIONS -------------------------------
t = t(:);

%% ==================== TEST 1: exact for polynomials =====================
fprintf('\n--- Test 1: exactness on polynomials ---\n');
fprintf('Lobatto nodes should be exact in the node-to-node integration up to nnodes-1 = %d\n', nnodes-1);
maxDeg  = 8;
poly_ok = true;
for k = 0:maxDeg
    f       = t.^k;
    Iexact  = (t(2:end).^(k+1) - t(1:end-1).^(k+1)) / (k+1);
    Iapprox = S*f;
    err     = max(abs(Iapprox - Iexact));
    pass    = err < tol_poly * max(1, abs(Iexact(end)));
    fprintf('  degree %2d : max abs error = %.3e  [%s]\n', k, err, tern(pass,'PASS','FAIL'));
    poly_ok = poly_ok && pass;
end

%% ==================== TEST 1b: exact for polynomials =====================
fprintf('\n--- Test 1b: exactness on polynomials ---\n');
fprintf('Lobatto nodes should be exact up to degree 2*nnodes-3 = %d\n', 2*nnodes-3);
maxDeg  = 8;
poly_ok = true;
for k = 0:maxDeg
    f       = t.^k;
    Iexact  = (t(end).^(k+1) - t(1).^(k+1)) / (k+1);
    Iapprox = sum(S*f);
    err     = max(abs(Iapprox - Iexact));
    pass    = err < tol_poly * max(1, abs(Iexact(end)));
    fprintf('  degree %2d : max abs error = %.3e  [%s]\n', k, err, tern(pass,'PASS','FAIL'));
    poly_ok = poly_ok && pass;
end

%% ==================== TEST 2: constant function ==========================
fprintf('\n--- Test 2: f(t) = 1  =>  S*f should equal (t - t0) ---\n');
f1  = ones(N,1);
err = max(abs(S*f1 - (t(2:end) - t(1:end-1))));
fprintf('  max abs error = %.3e  [%s]\n', err, tern(err < tol_poly, 'PASS', 'FAIL'));

%% ==================== TEST 3: smooth transcendental functions ===========
fprintf('\n--- Test 3: smooth functions vs high-accuracy quadrature ---\n');
mid       = (t0+t1)/2;
testFcns  = {@(x) sin(x), @(x) exp(0.3*x), @(x) 1./(1 + (x-mid).^2)};
names     = {'sin(t)', 'exp(0.3t)', '1/(1+(t-mid)^2)'};
for i = 1:numel(testFcns)
    fh = testFcns{i};
    F  = fh(t);
    Iapprox = S*F;
    Iexact  = zeros(N-1,1);
    for k = 2:N
        Iexact(k-1) = integral(fh, t(k-1), t(k), 'AbsTol', 1e-12, 'RelTol', 1e-12);
    end
    err = max(abs(Iapprox - Iexact));
    fprintf('  %-16s max abs error = %.3e  [%s]\n', names{i}, err, tern(err < tol_smooth, 'PASS', 'FAIL'));
end

%% ==================== TEST 4: convergence study ===========================
fprintf('\n--- Test 4: convergence as resolution increases ---\n');
fh          = @(x) exp(sin(x));
refineList  = 1:5;   % multiplies `intervals` by this factor -- edit as needed
errs        = zeros(size(refineList));
for j = 1:numel(refineList)
    ivals = intervals*refineList(j);
    ss    = (t1-t0)/ivals;
    [Sj, tj] = computeSpecMat(t0, t1, ss, ivals, nnodes, nodeType);
    tj      = tj(:);
    Fj      = fh(tj);
    Ij      = Sj*Fj;
    Iex     = arrayfun(@(a,b) integral(fh, a, b, 'AbsTol', 1e-13, 'RelTol', 1e-13), tj(1:end-1), tj(2:end));
    errs(j) = max(abs(Ij - Iex));
    fprintf('  intervals = %3d : max abs error = %.3e\n', ivals, errs(j));
end
% figure;
% semilogy(intervals*refineList, errs, '-o'); grid on;
% xlabel('number of sub-intervals'); ylabel('max abs error');
% title('Convergence of spectral integration matrix');

%%==================== TEST 5: telescoping / total-integral check =========
fprintf('\n--- Test 5: sum(S*f) should equal the integral over the whole domain ---\n');
for i = 1:numel(testFcns)
    fh = testFcns{i};
    F  = fh(t);
    total_approx = sum(S*F);
    total_exact  = integral(fh, t0, t1, 'AbsTol', 1e-12, 'RelTol', 1e-12);
    err = abs(total_approx - total_exact);
    fprintf('  %-16s total error = %.3e  [%s]\n', names{i}, err, tern(err < tol_smooth, 'PASS', 'FAIL'));
end
fprintf('  (cumsum(S*f) should also match a cumulative-from-t0 integration matrix,\n');
fprintf('   if you have one to compare against.)\n');
 


%% ==================== SUMMARY ==============================================
fprintf('\n=== SUMMARY ===\n');
fprintf('Polynomial exactness test: %s\n', tern(poly_ok, 'PASS', 'see above'));
fprintf('Review Test 3/4/5 above for smooth-function accuracy and convergence behavior.\n');

%% ------------------------- helper functions ---------------------------------
function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end

function t = getNodeLocations(t0, step_size, intervals, nodes, sharedBoundary)
    % Reconstructs physical collocation-node locations, assuming `nodes`
    % is either a count (defaulting to Legendre-Gauss-Lobatto points) or
    % an explicit vector of reference nodes on [-1,1], affinely mapped
    % onto each of the `intervals` sub-intervals of length step_size,
    % concatenated in order. EDIT to mirror your own implementation
    % exactly if it differs.
    if isnumeric(nodes) && isscalar(nodes)
        ref = lglnodes(nodes);      % n collocation points on [-1,1]
    elseif isnumeric(nodes) && isvector(nodes)
        ref = nodes(:);             % already the reference nodes
    else
        error('Edit getNodeLocations to handle nodes of class %s', class(nodes));
    end
    ref = sort(ref);
    t = [];
    for k = 0:intervals-1
        a  = t0 + k*step_size;
        b  = a + step_size;
        tk = a + (b-a)*(ref+1)/2;   % map [-1,1] -> [a,b]
        if sharedBoundary && k > 0
            tk(1) = [];             % drop duplicate shared left endpoint
        end
        t = [t; tk]; %#ok<AGROW>
    end
end

function x = lglnodes(n)
    % Legendre-Gauss-Lobatto nodes on [-1,1] (n points), via Newton
    % iteration on the Legendre polynomial recurrence.
    if n == 1
        x = 0; return;
    end
    N1 = n-1;
    x  = cos(pi*(0:N1)/N1)';    % Chebyshev-Gauss-Lobatto initial guess
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