% by Michael's python code
%Nbridge=100000000
%Nsteps=10000
%Bigdt = 6.25e-2
%meanRef =  0.5150808053464835; %0.5151244091857344; %0.5150843616737727;   %0.5150411119107488;  %second val: with 1000 Nsteps
%varRef  =   0.016527713576287767; %0.01653093956051733; %0.016524927225187978; %0.016534166369012896;
%latest
meanRef =  0.5151171113731648;
varRef =   0.01653079756143411;



% ODE45 for SBB-ODE, 
%1e5 samples
%meanODE(1,:) = [ 0.515034983898898,0.514983278911279, 0.514961723568478, 0.514933694338022 ];
%varODE(1,:)  = [0.016651891788958, 0.016648060400191, 0.016645849042841, 0.016643950340426];
%1e7 samples
% meanODE(1,:) = [ 0.515212563662107, 0.515163221025253, 0.515140914259654, 0.515112775339944,0.510261832621566];
% varODE(1,:)  = [ 0.016548384255890, 0.016545538145541, 0.016543375650266, 0.016541493358081,0.016657682581378];
% 
% errMean(1,:) = abs(meanODE(1,:) - meanRef) ./ meanRef; %[1.130759e-04, 6.600614e-05, 3.028814e-05, 7.423824e-05 ];
% errVar(1,:)  = abs(varODE(1,:)  - varRef)  ./ varRef;  %[2.283605e-04, 2.260874e-04, 3.126635e-04, 3.688699e-04 ];

% %DT=1/16; 5 col pts; m=1,2,5,10; 1e7 realizations
meanSDC(1,:) = [5.152126e-01, 5.151636e-01, 5.151936e-01, 5.152095e-01 ];
varSDC(1,:)  = [1.654838e-02, 1.654600e-02, 1.663490e-02, 1.667822e-02];
errMean(1,:) = abs(meanSDC(1,:) - meanRef) ./ meanRef; %[1.130759e-04, 6.600614e-05, 3.028814e-05, 7.423824e-05 ];
errVar(1,:)  = abs(varSDC(1,:)  - varRef)  ./ varRef;  %[2.283605e-04, 2.260874e-04, 3.126635e-04, 3.688699e-04 ];

% %DT=1/16; 7 col pts; m=1,2,5,10; 1e7 realizations
meanSDC(2,:) = [5.152126e-01];
varSDC(2,:)  = [1.654838e-02];
errMean(2,:) = abs(meanSDC(2,:) - meanRef) ./ meanRef; %[1.130759e-04, 6.600614e-05, 3.028814e-05, 7.423824e-05 ];
errVar(2,:)  = abs(varSDC(2,:)  - varRef)  ./ varRef;  %[2.283605e-04, 2.260874e-04, 3.126635e-04, 3.688699e-04 ];


% %DT=1/16; 5 col pts; m=1,2,5,10; 1e5 realizations
% meanSDC(1,:) = [5.150350e-01, 5.149836e-01, 5.150085e-01, 5.150369e-01];
% varSDC(1,:)  = [1.665189e-02, 1.664851e-02, 1.673967e-02, 1.678588e-02];
% errMean(1,:) = abs(meanSDC(1,:) - meanRef) ./ meanRef; %[1.130759e-04, 6.600614e-05, 3.028814e-05, 7.423824e-05 ];
% errVar(1,:)  = abs(varSDC(1,:)  - varRef)  ./ varRef;  %[2.283605e-04, 2.260874e-04, 3.126635e-04, 3.688699e-04 ];
% 
% %DT=1/16; 9 col pts; m=1,2,5,10; 1e5 realizations
% meanSDC(2,:) = [5.150350e-01, 5.149833e-01, 5.149616e-01, 5.149736e-01 ];
% varSDC(2,:)  = [1.665189e-02, 1.664806e-02, 1.664622e-02, 1.669423e-02 ];
% errMean(2,:) = abs(meanSDC(2,:) - meanRef) ./ meanRef; 
% errVar(2,:)  = abs(varSDC(2,:)  - varRef)  ./ varRef; 
% 
% %DT=1/16; 15 col pts; m=1,2,5,10; 1e5 realizations
% meanSDC(3,:) = [ 5.150350e-01, 5.149833e-01, 5.149619e-01, 5.149613e-01  ];
% varSDC(3,:) = [1.665189e-02,1.664806e-02,1.664587e-02,1.665228e-02   ];
% errMean(3,:) = abs(meanSDC(3,:) - meanRef) ./ meanRef; 
% errVar(3,:)  = abs(varSDC(3,:)  - varRef)  ./ varRef; 
% %errMean(3,:) = [5.150350e-01 , 5.149833e-01, 5.149619e-01  ];
% %errVar(3,:)  = [1..665189e-02 , 1.664806e-02 , 1.664587e-02];

m=[1,2,5,10];

figure();
plot(m, errMean(1,:), '-o', 'LineWidth',2, 'MarkerSize',5);
% legend('Matlab ODE45');
xticks(m);
xlabel('# bridge terms');
ylabel('relative error in mean');
title('Nonlin example, 1 timestep,  1e7 samples');
saveas(gcf, 'TP3-SDC-relerr-mean.eps', 'epsc');

figure();
plot(m, errVar(1,:), '-o', 'LineWidth',2, 'MarkerSize',5);
% legend('Matlab ODE45');
xticks(m);
xlabel('# bridge terms');
ylabel('relative error in variance');
title('Nonlin example, 1 timestep, 1e7 samples');
saveas(gcf, 'TP3-SDC-relerr-var.eps', 'epsc');


%m, errMean(3,:), '-s',
% figure();
% plot(m, errMean(1,:), '-o', m, errMean(2,:), '-x', m, errMean(3,:), '-s','LineWidth',2, 'MarkerSize',5);
% legend('5 Lobatto nodes', '9 Lobatto nodes', '15 Lobatto nodes');
% xticks(m);
% xlabel('# bridge terms');
% ylabel('relative error in mean');
% title('Nonlinear example, 1 timestep only,  1e5 samples');
% saveas(gcf, 'TP3-relerr-mean.eps', 'epsc');
% 
% figure();
% plot(m, errVar(1,:), '-o', m, errVar(2,:), '-x', m, errVar(3,:), '-s', 'LineWidth',2, 'MarkerSize',5);
% legend('5 Lobatto nodes', '9 Lobatto nodes', '15 Lobatto nodes');
% xticks(m);
% xlabel('# bridge terms');
% ylabel('relative error in variance');
% title('Nonlinear example, 1 timestep only, 1e5 samples');
% saveas(gcf, 'TP3-relerr-var.eps', 'epsc');
