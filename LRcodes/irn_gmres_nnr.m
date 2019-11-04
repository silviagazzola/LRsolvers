function [Xfinals, Xbest, Enrm_tot, Rnrm_tot, cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_gmres_nnr(A, b, xex, x0, parameters)


cycles = parameters.cycles;
thr = parameters.thr;
reg = parameters.reg;
p = parameters.p;
maxIt = parameters.maxIt;
eta = parameters.eta;
nl = parameters.nl;
weigthtype = parameters.weigthtype;
thrstop = parameters.thrstop;
% weigthtype = 'trunc'; % weigthtype = 'sqrt'; % 

N = length(x0);
n = sqrt(N);
totIt = cycles*maxIt;
Xfinals = zeros(N,cycles);
Xbest = zeros(N,cycles);
Enrm_tot = zeros(totIt, 1);
Rnrm_tot = zeros(totIt, 1);
Lambda_tot = zeros(totIt, 1);
cyclesIt = zeros(cycles, 1);
sXcycles = zeros(n, cycles);
cyclesBest = zeros(cycles, 1);

% first cycle of iterations (with no regularization)
optnnr.p = p;
optnnr.maxIt = maxIt;
optnnr.reg = reg;
optnnr.eta = eta;
optnnr.nl = nl;

UX = eye(n); VX = eye(n); SX = ones(n); 
totIt_temp = 0;

for i = 1:cycles
    optnnr.svd.U = UX; optnnr.svd.V = VX; optnnr.svd.S = SX;
    [X_nnr, Enrm_nnr, Rnrm_nnr, Lambda_nnr] = gmres_nnr(A, b, xex, x0, optnnr);
    [~, indmin] = min(Enrm_nnr(Enrm_nnr ~= 0));
    % trying to set a stopping criterion (after all the iterations are performed)
    Rnrm_nnr = Rnrm_nnr(Rnrm_nnr ~= 0);
    ind = find(Rnrm_nnr <= eta*nl, 1);
    if isempty(ind)
        ind = find(Enrm_nnr == 0, 1);
        ind = ind - 1;
    end
    if isempty(ind)
        ind = length(Rnrm_nnr);
    end
    if isstruct(reg)
        reg.value = Lambda_nnr(ind);
        optnnr.reg = reg;
    end
    Xtemp = X_nnr(:, ind);
    Xfinals(:,i) = Xtemp;
    Xbest(:,i) = X_nnr(:, indmin);
    cyclesBest(i) = indmin;
    cyclesIt(:,i) = ind;
    Enrm_tot(totIt_temp+1:totIt_temp+ind) = Enrm_nnr(1:ind);
    Rnrm_tot(totIt_temp+1:totIt_temp+ind) = Rnrm_nnr(1:ind);
    Lambda_tot(totIt_temp+1:totIt_temp+ind) = Lambda_nnr(1:ind);
    % setting the next regularization matrix
    [UX, SX, VX] = svd(reshape(Xtemp, n, n));
    sX = diag(SX);
    sXcycles(:,i) = sX./sX(1);
    if strcmp(weigthtype, 'trunc')
        diffsX = sX(1:end-1)-sX(2:end);
        diffsXr = diffsX./diffsX(1);
        % diffsXr = diffsX./sX(1:end-1);
        trunc = find(diffsXr<thr, 1);
        %     if trunc > 1
        %         trunc = trunc - 1;
        %     end
        SX(trunc:end,trunc:end) = 0;
    end
    if i>2 && norm(sXcycles(:,i-1)-sXcycles(:,i))<thrstop
        break
    end
    totIt_temp = totIt_temp + ind;
end

totIt = totIt_temp;

Enrm_tot = Enrm_tot(1:totIt);
Rnrm_tot = Rnrm_tot(1:totIt);

    
