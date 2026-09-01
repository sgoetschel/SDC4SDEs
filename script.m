close all
clear all
% %%problem definition
% %d = 'const';
% %d = 'exp';
% %d = 'TP2';
% %d = 'TP3';
% %d = 'TP4';
% 
% %%method choice
% %sdc_choice = 'SDC_BB';
% %sdc_choice = 'EM';
% %sdc_choice = 'Mil';
% %sdc_choice = 'implMilDampNM';
% %sdc_choice = 'implEMDampNM';
% 
% collocation nodes type
% %nodes ='equidist';
% %nodes ='legendre';
% %nodes ='lobatto';
% %includes right border
% %nodes ='radau';
% 
% problem initialization
% %strInit = 'constInit';
% %strInit = 'eulerSBBInit';
% %strInit = 'eMItoSDE';
% %strInit = 'eMItoSDE_implDampNM; %-->WORKS ONLY FOR TP 2!
%
% %different  main functions
% main_SDCSDEs_SweepsOrder
% main_SDCSDEs
% main_realIter_moments
%
% %%choose number of collocation points on sub time interval
% %colpoints = [3, 4, 5];
% 
% %max SDC sweeps
% %maxIter = [4, 6, 8];
% 
% %%variables
% %number of realizations
% %realIter = 1e3;
% 
% %Brownian bridge - m ~ fourier series terms
% %mBB = 1;%[1, 2, 3, 4, 5];
% 
% %%output information
% %plot_Sol = false;
% %plot_error = true;
%
% %exact/ref solution using SBB-ODE
% %SBB = true;
% 
% %starts with 2 intervals and has to be a multiple of 2
% %steps = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625];%, 0.0078125, 0.0078125/2];
% 
% % variables:(probDef, sde_solver, nodes,  strInit,  colpoints, maxIter, steps, NNfinest, realIter, m, plot_Sol, plot_error, SBB_ODE)            


%% const 
% problem const is not working properly!!!

%% TP 3 - 4 different initializations, and EM and Milstein
%fprintf(1, '\n problem: TP3 with lambda=1, beta=1 using EulerSBBInit; SDC with eulerSBBInit \n');
% deltaT = [1.0/32, 1.0/64, 1.0/128, 1.0/256, 1.0/512];
% %deltaT = [1.0/256];
% for m= [1,2,3,4,5]
%   for n=1:3
%     main_SDCSDEs('TP3',    'SDC_BB',   'lobatto', 'eulerSBBInit', 2+n, 100, deltaT, 512, 1e4,  m, false,   false, false);
%   end
% end

warning('off','MATLAB:MKDIR:DirectoryExists');
diary myLogFile
nSamples = 1e4;
for deltaT = [1/16; 1/32; 1/64; 1/128; 1/256; 1/512; 1/1024]
  for n=[4] % collocation points
    for m=[1] % Brownian Bridge terms
      %fprintf("\n\n=====================================================================================================\n");
      %fprintf(" n = %d, m = %d\n", n, m);
      
      %           (    d,  sde_solver,      nodes,       strInit, colpoints, maxIter,  steps, NNfinest, realIter, mBB, plot_Sol, plot_Error, SBB)
      main_SDCSDEs('TP2',    'SDC_BB',   'lobatto', 'eulerSBBInit',        n,      10, deltaT,     1024, nSamples,   m,    false,      true, false);
      %fprintf("\n\n=====================================================================================================\n");
    end
  end
end
diary off

% deltaT = [2^(-1), 2^(-2)]%, 2^(-3), 2^(-4), 2^(-5), 2^(-6)]; % time interval [0,2] -> SET NNFINEST CORRECTLY! (double it compared to other exampels)
% %m=2;
% n=3;
% nSamples=1e5;
% for m= [1,2,3,4,5]
%   for n=[1,2,3,4,5]
%     main_SDCSDEs('Mattingly',    'SDC_BB',   'lobatto', 'eulerSBBInit', 2+n, 100, deltaT, 128, nSamples,  m, false,   false, false);
%   end
% end


% deltaT = [1.0/256];
% for m= [1,2,3,4,5]
%   for n=1:7
%     main_SDCSDEs('exp',    'SDC_BB',   'lobatto', 'eulerSBBInit', 2+n, 100, deltaT, 4096, 1e5,  m, false,   false, true);
%   end
% end


%%% TP1 - 3 different initializations, and EM and Milstein 
%for n=n=1:3
% main_SDCSDEs('exp',    'SDC_BB',   'lobatto', 'constInit', 2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e3,       1, false,    false, false);
%
%  main_SDCSDEs('exp',    'SDC_BB',   'lobatto', 'eulerSBBInit',  2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e3,       1, false,  false, false);
% 
% main_SDCSDEs('exp',    'SDC_BB',   'lobatto', 'eMItoSDE', 2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e3,       1, false,    false, false);
%end
 
  

% solver = ["EM", "Mil"];
% parfor k=1:length(solver)
% main_SDCSDEs('exp',    solver(k),   'equidist', 'constInit', 3, 4, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e4,       1, false,    false, false);
% end
% 
% %% TP 2 - 3 different initializations, and EM and Milstein
% parfor n=1:3
% main_SDCSDEs('TP2',    'SDC_BB',   'lobatto', 'constInit', 2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% 
% main_SDCSDEs('TP2',    'SDC_BB',   'lobatto', 'eulerSBBInit', 2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% 
% main_SDCSDEs('TP2',    'SDC_BB',   'lobatto', 'eMItoSDE_implDampNM', 2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% main_SDCSDEs('TP2',    'SDC_BB',   'lobatto', 'eMItoSDE', 2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% 
% end
% 
% solver = ["implEMDampNM", "implMilDampNM"];
% parfor k=1:length(solver)
% main_SDCSDEs('TP2',    solver(k),   'equidist', 'constInit', 3, 4, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% end
% 
% %% TP 4 - 3 different initializations, and EM and Milstein
% parfor n=1:3
% main_SDCSDEs('TP4',    'SDC_BB',   'lobatto', 'constInit',  2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% 
% main_SDCSDEs('TP4',    'SDC_BB',   'lobatto', 'eulerSBBInit',  2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% 
% main_SDCSDEs('TP4',    'SDC_BB',   'lobatto', 'eMItoSDE',  2+n, 2+2*n, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% end 
% 
% solver = ["EM", "Mil"];
% parfor k=1:length(solver)
% main_SDCSDEs('TP4',    solver(k),   'equidist', 'constInit', 3, 4, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
%     1e5,       1, false,    false, false);
% %main_SDCSDEs('TP4',    'Mil',   'equidist', 'constInit', 3, 4, [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625e-04, 0.0078125/16, 0.0078125/32, 0.0078125/64, 0.0078125/128], 16384, ...
% %    1e5,       1, false,    false, false);
% end


%%% Test for sweeps vs order - calls main_SDCSDEs_SweepsOrder
%for n=1:9
%main_SDCSDEs_SweepsOrder('exp', 'SDC_BB', 'lobatto', 'eMItoS', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);
%main_SDCSDEs_SweepsOrder('exp', 'SDC_BB', 'lobatto', 'eulerSBBInit', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);

%main_SDCSDEs_SweepsOrder('TP2', 'SDC_BB', 'lobatto', 'eulerSBBInit', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);
%main_SDCSDEs_SweepsOrder('TP2', 'SDC_BB', 'lobatto', 'implDampNM_eMItoSDE', 5, n-1 , [2^-3, 2^-4, 2^-5, 2^-6], 1e4, 1, false, false, false);

%main_SDCSDEs_SweepsOrder('TP3', 'SDC_BB', 'lobatto', 'eulerSBBInit', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);

%main_SDCSDEs_SweepsOrder('TP4', 'SDC_BB', 'lobatto', 'eulerSBBInit',5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);
%main_SDCSDEs_SweepsOrder('TP4', 'SDC_BB', 'lobatto', 'eMItoS', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);

%end

%parfor n=1:10
%main_SDCSDEs_SweepsOrder('exp', 'SDC_BB', 'lobatto', 'constInit', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);

%main_SDCSDEs_SweepsOrder('TP2', 'SDC_BB', 'lobatto', 'constInit', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);

%main_SDCSDEs_SweepsOrder('TP4', 'SDC_BB', 'lobatto', 'constInit', 5, n-1, [2^-3, 2^-4, 2^-5, 2^-6], 1024, 1e4, 1, false, false, false);

%end
