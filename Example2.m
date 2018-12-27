%% Example2.m
%
% This code illustrates how to include a transformation matrix in the codes 
% corresponding to the paper
%     "Flexible Krylov Methods for l_p Regularization"
%       - Chung and Gazzola, 2019

IRtools_setup
rng(0)
clear, clc
%% create data
x_true = double(imresize(imread('cameraman.tif'),.5));
n = size(x_true,1);

% define the PSF
x1 = -fix(n/2):ceil(n/2)-1;
[X1,Y1] = meshgrid(x1,x1);
PSF = (X1.^2 + Y1.^2) <= 4^2;
PSF = PSF/sum(PSF(:));

optblur.trueImage = reshape(x_true, n, n);
optblur.PSF = PSF;
optblur.BC = 'reflective';
optblur.CommitCrime = 'on';
[A, b, x_true, ProbInfo] = PRblur(optblur);
    
% Use wavelet transform
wname = 'db1';
wlevels = 3;
Trans = FreqMatrix('dwt', [n n], wname, wlevels);
w = Trans*x_true;
fprintf('Image domain: %d \n',sum(x_true(:)<eps))
fprintf('Wavel domain: %d \n',sum(abs(w(:))<eps))

nl = .01; % Noise level
bn = PRnoise(b, nl);

figure, subplot(1,3,1), imshow(reshape(x_true,n,n),[]), title('true')
im = viswavdec(w,Trans,2); im(abs(im)<eps) = 0;
subplot(1,3,2), imshow(im,[]), title('true (wavelet)')
subplot(1,3,3), imshow(reshape(bn,n,n),[]), title('observed')

Afcn = @(xx,tflag) OPblur_transform(xx, A, Trans, tflag);
bntransf = Trans*bn;

opt.SparsityTrans = 'dwt';
opt.wlevels = 3;
opt.x_true = x_true;
opt.NoStop = 'on';
opt.RegParam = 'off';
opt.NoiseLevel = nl;
maxit=150;
K = [1, 10:10:maxit];
%% run LSQR and FLSQR (no regularization)
[X_LSQR, info_LSQR] = IRhybrid_lsqr(A, bn, K, opt);
[X_FLSQR, info_FLSQR] = IRhybrid_flsqr(Afcn, bntransf, K, opt);
%% run FLSQR-I
opt = IRset(opt, 'RegParam', 'discrep');
[X_FLSQRi, info_FLSQRi] = IRhybrid_flsqr(Afcn, bntransf, K, opt);
%% run FLSQR-R
opt = IRset(opt, 'hybridvariant', 'R');
% secant update parameter choice
[X_FLSQRr, info_FLSQRr] = IRhybrid_flsqr(Afcn, bntransf, K, opt);
%% run FAT-I
opt = IRset(opt, 'hybridvariant', 'I');
[X_FAT, info_FAT] = IRhybrid_fgmres(Afcn, bntransf, K, opt);
%% run FISTA
RegP = (info_FLSQRr.StopReg.RegP)^2;
opt.RegParam = RegP;
opt.shrink = 'on';
opt.xMin = -Inf;
[X_FISTA, info_FISTA] = IRfista(A, bn, K, opt);

%% Display images and plots
figure, lw = 2;
plot(info_FLSQR.Enrm, '-k','LineWidth',lw), hold on
plot(info_FLSQRi.Enrm, 'b-.','LineWidth',lw)
plot(info_FLSQRr.Enrm, 'r--','LineWidth',lw)
plot(info_LSQR.Enrm, '-*c','LineWidth',2,'MarkerIndices',1:10:length(info_LSQR.Enrm),'MarkerSize',12), hold on
plot(info_FLSQRi.StopReg.It, info_FLSQRi.Enrm(info_FLSQRi.StopReg.It),'b*', 'MarkerSize',16, 'LineWidth',lw)
plot(info_FLSQRr.StopReg.It, info_FLSQRr.Enrm(info_FLSQRr.StopReg.It),'rd', 'MarkerSize',16, 'LineWidth',lw)
xlabel('Iteration'), ylabel('Relative Error')
legend('FLSQR', 'FLSQR-I','FLSQR-R','LSQR')
axis([0,maxit,.005, .4])

% Compare to other methods
figure, 
semilogy(info_FLSQRr.Enrm, '--r','LineWidth',lw), hold on
semilogy(info_FAT.Enrm, '-.m','LineWidth',lw)
semilogy(info_FISTA.Enrm, '-oc','LineWidth',2,'MarkerIndices',1:10:length(info_FISTA.Enrm),'MarkerSize',12)
semilogy(info_FLSQRr.StopReg.It, info_FLSQRr.Enrm(info_FLSQRr.StopReg.It),'rd', 'MarkerSize',16, 'LineWidth',lw)
semilogy(info_FAT.StopReg.It, info_FAT.Enrm(info_FAT.StopReg.It),'ms', 'MarkerSize',16, 'LineWidth',lw)
xlabel('Iteration'), ylabel('Relative Error')
legend('FLSQR-R', 'GAT', 'FISTA')
axis([0,maxit,.07,.3])
%% Show images
xbest_lsqr = info_LSQR.BestReg.X;
xbest_flsqr = info_FLSQR.BestReg.X;
xbest_flsqrr = info_FLSQRr.BestReg.X;
xbest_flsqri = info_FLSQRi.BestReg.X;
Xbest_lsqr = reshape(xbest_lsqr, n, n);
Xbestt_lsqr = Xbest_lsqr(33:33+127, 30:30+127);
Xbest_flsqr = reshape(xbest_flsqr, n, n);
Xbestt_flsqr = Xbest_flsqr(33:33+127, 30:30+127);
Xbest_flsqrr = reshape(xbest_flsqrr, n, n);
Xbestt_flsqrr = Xbest_flsqrr(33:33+127, 30:30+127);
Xbest_flsqri = reshape(xbest_flsqri, n, n);
Xbestt_flsqri = Xbest_flsqri(33:33+127, 30:30+127);

x_true = reshape(x_true, n, n);
e_lsqr = abs(Xbest_lsqr - x_true);
et_lsqr = e_lsqr(33:33+127, 30:30+127);
e_flsqr = abs(Xbest_flsqr - x_true);
et_flsqr = e_flsqr(33:33+127, 30:30+127);
e_flsqrr = abs(Xbest_flsqrr - x_true);
et_flsqrr = e_flsqrr(33:33+127, 30:30+127);
e_flsqri = abs(Xbest_flsqri - x_true);
et_flsqri = e_flsqri(33:33+127, 30:30+127);

figure, 
subplot(2,4,1), imagesc(Xbestt_lsqr), colormap gray, axis square, axis image, axis off, title('LSQR')
subplot(2,4,2), imagesc(Xbestt_flsqr), colormap gray, axis square, axis image, axis off, title('FLSQR')
subplot(2,4,3), imagesc(Xbestt_flsqri), colormap gray, axis square, axis image, axis off, title('FLSQR-I')
subplot(2,4,4), imagesc(Xbestt_flsqrr), colormap gray, axis square, axis image, axis off, title('FLSQR-R')
subplot(2,4,5), imagesc(1-et_lsqr/max(max(et_lsqr))), colormap gray, axis square, axis image, axis off
subplot(2,4,6), imagesc(1-et_flsqr/max(max(et_flsqr))), colormap gray, axis square, axis image, axis off
subplot(2,4,7), imagesc(1-et_flsqri/max(max(et_flsqri))), colormap gray, axis square, axis image, axis off
subplot(2,4,8), imagesc(1-et_flsqrr/max(max(et_flsqrr))), colormap gray, axis square, axis image, axis off
