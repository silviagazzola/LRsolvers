%% Example1.m
%
% This code illustrates how to use the codes corresponding to the paper
%     "Flexible Krylov Methods for l_p Regularization"
%       - Chung and Gazzola, 2019

IRtools_setup
rng(0)
clear, clc
%% create data
n = 256;
load star_testpb

optblur.trueImage = reshape(x_true, n, n);
optblur.PSF = PSF;
optblur.BC = 'zero';
optblur.CommitCrime = 'on';
[A, b, x_true, ProbInfo] = PRblur(optblur);
W = 1; % No transformation included
nl = 5e-2; % Noise level
bn = PRnoise(b, nl); % Observed Image

options = IRset('x_true', x_true(:), 'NoStop', 'on', 'RegParam', 0,...
    'NoiseLevel', nl, 'DecompOut', 'on');
maxit = 150;
K = 1:maxit;
%% run FLSQR and LSQR (no regularization)
[X_FLSQR, info_FLSQR] = IRhybrid_flsqr(A, bn, K, options);
[X_LSQR, info_LSQR] = IRhybrid_lsqr(A, bn, K, options);
%% run FLSQR-I with different parameter choice methods
% secant update parameter choice
options = IRset(options, 'RegParam', 'discrep');
[X_FLSQRi, info_FLSQRi] = IRhybrid_flsqr(A, bn, K, options);
% classical discrepancy principle
options = IRset(options, 'RegParam', 'discrepit');
[X_FLSQRidp, info_FLSQRidp] = IRhybrid_flsqr(A, bn, K, options);
%% run FLSQR-R with different parameter choice methods
options = IRset(options, 'RegParam', 'discrep', 'hybridvariant', 'R');
% secant update parameter choice
[X_FLSQRr, info_FLSQRr] = IRhybrid_flsqr(A, bn, K, options);
% classical discrepancy principle
options = IRset(options, 'RegParam', 'discrepit');
[X_FLSQRrdp, info_FLSQRrdp] = IRhybrid_flsqr(A, bn, K, options);
% optimal regularization parameter
options = IRset(options, 'RegParam', 'optimal');
[X_FLSQRropt, info_FLSQRropt] = IRhybrid_flsqr(A, bn, K, options);
%% compare with other solvers
RegP = info_FLSQRr.StopReg.RegP;
% run GAT
options.RegParam = 'discrep';
[X_FAT, info_FAT] = IRhybrid_fgmres(A, bn, K, options);
% run FISTA
options.shrink = 'on';
options.RegParam = RegP^2;
options.xMin = -Inf;
[X_FISTA, info_FISTA] = IRfista(A, bn, K, options);

%% Display Images and Plots
figure, subplot(1,3,1), imshow(reshape(x_true,n,n),[]), title('true')
subplot(1,3,2), imshow(reshape(PSF,n,n),[]), title('PSF')
subplot(1,3,3), imshow(reshape(bn,n,n),[]), title('observed')

figure, lw = 2; 
plot(info_FLSQR.Enrm, '-k','LineWidth',lw), hold on
plot(info_FLSQRi.Enrm, '-.b','LineWidth',lw)
plot(info_FLSQRr.Enrm, '--r','LineWidth',lw)
plot(info_LSQR.Enrm, '-*c','LineWidth',2,'MarkerIndices',1:10:length(info_LSQR.Enrm),'MarkerSize',12)
plot(info_FLSQRi.StopReg.It, info_FLSQRi.Enrm(info_FLSQRi.StopReg.It),'b*', 'MarkerSize',16, 'LineWidth',lw)
plot(info_FLSQRr.StopReg.It, info_FLSQRr.Enrm(info_FLSQRr.StopReg.It),'rd', 'MarkerSize',16, 'LineWidth',lw)
xlabel('Iteration'), ylabel('Relative Error')
legend('FLSQR', 'FLSQR-I','FLSQR-R','LSQR')
axis([0,maxit,0,1])

figure,
plot(info_FLSQR.Enrm, '-k','LineWidth',lw), hold on
plot(info_FLSQRr.Enrm, '--r','LineWidth',lw)
plot(info_FLSQRrdp.Enrm, ':om','LineWidth',lw,'MarkerIndices',1:10:maxit,'MarkerSize',10)
plot(info_FLSQRropt.Enrm, '-.*b','LineWidth',lw,'MarkerIndices',1:10:maxit,'MarkerSize',10)
xlabel('Iteration'), ylabel('Relative Error')
legend('FLSQR', 'FLSQR-R', 'FLSQR-R dp', 'FLSQR-R opt')
axis([0,maxit,0,1])

figure,
plot(info_FLSQRr.Enrm, '--r','LineWidth',lw), hold on
plot(info_FAT.Enrm, '-.m','LineWidth',lw)
plot(info_FISTA.Enrm, '-oc','LineWidth',2,'MarkerIndices',1:10:length(info_FISTA.Enrm),'MarkerSize',12)
plot(info_FLSQRr.StopReg.It, info_FLSQRr.Enrm(info_FLSQRr.StopReg.It),'rd', 'MarkerSize',16, 'LineWidth',lw)
plot(info_FAT.StopReg.It, info_FAT.Enrm(info_FAT.StopReg.It),'ms', 'MarkerSize',16, 'LineWidth',lw)
xlabel('Iteration'), ylabel('Relative Error')
legend('FLSQR-R', 'GAT', 'FISTA')
axis([0,maxit,0,1])

%% Show the basis vectors
Z = info_FLSQRrdp.Z; V = info_LSQR.V;
figure, 
subplot(2, 3, 1), imshow(reshape(Z(:,10), n, n), []), title('FLSQR-R, k = 10')
subplot(2, 3, 2), imshow(reshape(Z(:,20), n, n), []), title('FLSQR-R, k = 20')
subplot(2, 3, 3), imshow(reshape(Z(:,100), n, n), []), title('FLSQR-R, k = 100')
subplot(2, 3, 4), imshow(reshape(V(:,10), n, n), []), title('LSQR, k = 10')
subplot(2, 3, 5), imshow(reshape(V(:,20), n, n), []), title('LSQR, k = 20')
subplot(2, 3, 6), imshow(reshape(V(:,100), n, n), []), title('LSQR, k = 100')
