clear, clc

%% first generate the test problem
optbl = PRblur('defaults');
optbl = PRset(optbl, 'trueImage', 'dot2', 'BlurLevel', 'medium');
% Medium Gaussian blur
[A, b, x, ProbInfo] = PRblur(optbl);

N = length(x); n = sqrt(N);

nl = 1e-3; % 
rng = 0; % setting the random number generator seed -- to obtain reproducible results
bn = PRnoise(b, nl);

%% SOLVERS BASED ON THE ARNOLDI ALGORITHM

% standard GMRES (as implemented in IR Tools)
maxIt = 100;
opth.RegParam = 0; % purely iterative version (no hybrid)
opth.NoiseLevel = nl;
opth.NoStop = 'on';
opth.x_true = x;
opth.DecompOut = 'on';

[X_gmres, info_gmres] = IRhybrid_gmres(A, bn, 1:maxIt, opth);

% IRN-GMRES-NNR
eta = 1.5; % safety threshold for the discrepancy principle
x0 = zeros(N, 1);
p = 1;

parameters.cycles = 5; % outer iterations
parameters.thr = 1e-3;
parameters.reg = 0; % parameters.reg = 'discrep'; % 
parameters.p = p;
parameters.maxIt = maxIt; % inner iterations
parameters.eta = eta;
parameters.nl = nl;
parameters.weigthtype = 'sqrt'; % preconditionong type
parameters.thrstop = 1e-8; % for the stopping criterion for the outer iterations

[Xfinals_irn, Xbest_irn, Enrm_tot_irn, Rnrm_tot_irn, cyclesIt_irn, sXcycles_irn, cyclesBest_irn] = irn_gmres_nnr(A, bn, x, x0, parameters);
 
% FGMRES-NNR

optnnr.p = p;
optnnr.maxIt = 200;
optnnr.regmat = 'I';
optnnr.reg = 0; % purely iterative version (no hybrid)
optnnr.eta = eta;
optnnr.nl = nl;
optnnr.svdbasis = 1;
[X_fnnr,Enrm_fnnr,Rnrm_fnnr] = fgmres_nnr(A, bn, x, x0, optnnr);

% LR-FGMRES
[X_LRgm,Rnrm_LRgm,Enrm_LRgm] = FGMRES_LRP(A,bn,maxIt,x0,2,2,x,1e-4,1,0,1.1,nl);

% RS-LR_GMRES
[X_RS, RelRes_RS, RelErr_RS] = RS_GMRES_LRP(A,bn,10,20,x0,2,2,x);

%% Comparison with SVT
[X_SVT,RelErr_SVT,RelRes_SVT] = SVT(A,bn,x,1,100,1e-3,2);

%% Displaying the results

% History of the relative errors
figure, semilogy(info_gmres.Enrm, 'LineWidth', 2), hold on
semilogy(Enrm_tot_irn, 'LineWidth', 2)
semilogy(Enrm_fnnr, 'LineWidth', 2)
semilogy(Enrm_LRgm, 'LineWidth', 2)
semilogy(RelErr_RS, 'LineWidth', 2)
semilogy(RelErr_SVT, 'LineWidth', 2)
legend('GMRES', 'IRN-GMRES-NNR', 'FGMRES-NNR', 'LR-FGMRES', 'RS-LR-GMRES', 'SVT')
xlabel('(Total) Iteration Count')
ylabel('Relative Error')

% Best IRN-GMRES-NNR solution
Xb_irnbest = reshape(Xbest_irn(:,5), n, n);
% IRN-GMRES-NNR solution at the end of the 1st cycle
Xb_irn1 = reshape(Xfinals_irn(:,1), n, n); % 
% IRN-GMRES-NNR solution at the end of the 2nd cycle
Xb_irn2 = reshape(Xfinals_irn(:,2), n, n);
% best FGMRES-NNR solution
[~, ind] = min(Enrm_fnnr);
Xb_fnnr = reshape(X_fnnr(:,ind), n, n);

% Plotting the data
figure
subplot(1,2,1), imagesc(reshape(x,256,256)), axis image, axis off, title('exact')
subplot(1,2,2), imagesc(reshape(bn,256,256)), axis image, axis off, title('corrupted')

% Plotting some of the best solutions
figure
subplot(2,2,1), imagesc(Xb_irnbest), axis image, axis off, title('IRN-GMRES-NNR, best')
subplot(2,2,2), imagesc(Xb_fnnr), axis image, axis off, title('FGMRES-NNR, best')
subplot(2,2,3), imagesc(Xb_irn1), axis image, axis off, title('IRN-GMRES-NNR, 1st cycle')
subplot(2,2,4), imagesc(Xb_irn2), axis image, axis off, title('IRN-GMRES-NNR, 2nd cycle')

