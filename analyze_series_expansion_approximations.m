clear all, close all, clc

%----------------------------------------------------------------------
% System definition
%----------------------------------------------------------------------
% first second-order system
G_dc1 = 2;
fn1 = 50;
varsigma1 = 0.005;
omega_n1 = 2*pi*fn1;

b1 = -(1i)*G_dc1*omega_n1^1/(2*sqrt(1-varsigma1^2));
p1 = -varsigma1*omega_n1+(1i)*omega_n1*sqrt(1-varsigma1^2);

% another second-order system
G_dc2 = 2;
fn2 = 35;
varsigma2 = 0.005;
omega_n2 = 2*pi*fn2;

b2 = -(1i)*G_dc2*omega_n2^1/(2*sqrt(1-varsigma2^2));
p2 = -varsigma2*omega_n2+(1i)*omega_n2*sqrt(1-varsigma2^2);

% And yet another second-order system
G_dc3 = 2;
fn3 = 45;
varsigma3 = 0.005;
omega_n3 = 2*pi*fn3;

b3 = -(1i)*G_dc3*omega_n3^1/(2*sqrt(1-varsigma3^2));
p3 = -varsigma3*omega_n3+(1i)*omega_n3*sqrt(1-varsigma3^2);

% Normalize DC value
G1_dc = b1./(0-p1) + conj(b1)./(0-conj(p1));
G2_dc = b2./(0-p2) + conj(b2)./(0-conj(p2));
G3_dc = b3./(0-p3) + conj(b3)./(0-conj(p3));
G0_dc = G1_dc+G2_dc+G3_dc;
scale_dc = 1/abs(G0_dc);

%----------------------------------------------------------------------
% Relevant orders in series expansion approximation
%----------------------------------------------------------------------
R = 1;      % order of polynomial
nn = 2*R;   % Local bandwidth

%----------------------------------------------------------------------
% Data acquisition parameters 
%----------------------------------------------------------------------
N_mult = [1,2:0.5:10,10:10:100,100:100:1000];    % multiplier for data length in
                                                 % analysis of cut-off error
% Allocate matrices to store results
dsigma1 = zeros(1,length(N_mult));
dsigma2 = zeros(1,length(N_mult));
G_err1 = zeros(1,length(N_mult));
G_err2 = zeros(1,length(N_mult));

for ii=1:length(N_mult)
    %----------------------------------------------------------------------
    % Data acquisition parameters cont.
    %----------------------------------------------------------------------
    
    M = 16;                 % number of repeated experiments
    Np = 2^8;               % block length
    N = N_mult(ii)*M*Np;    % data length
    fs = Np;                % sampling frequency s.t. 1 period = 1 s
    df = fs/N;
    lines_N = [0:1:N-1]';
    freq_N = lines_N*df;
    Omega_N = (1i)*2*pi*freq_N;
    
    %----------------------------------------------------------------------
    % Compute transfer function G0
    %----------------------------------------------------------------------
    % Compute G0 on the current frequency grid
    G1 = b1./(Omega_N-p1) + conj(b1)./(Omega_N-conj(p1));
    G2 = b2./(Omega_N-p2) + conj(b2)./(Omega_N-conj(p2));
    G3 = b3./(Omega_N-p3) + conj(b3)./(Omega_N-conj(p3));
    G0_N = G1+G2+G3;
    G0_N = scale_dc*G0_N;
    
    %----------------------------------------------------------------------
    % Change of variables, substition to local variable
    %----------------------------------------------------------------------
    varsigma = (1i)*2*pi*[-nn:1:nn]*df;
    omega0 = 2*pi*fn1;
    
    p_sys = [p1 p2 p3];
    res_sys = scale_dc*[b1 b2 b3];
    alpha = real([p_sys, conj(p_sys)]);
    beta = imag([p_sys, conj(p_sys)]);
    pfe.residues = [res_sys, conj(res_sys)];
    pfe.poles = alpha+(1i)*(beta-omega0);
    
    %----------------------------------------------------------------------
    % Check conditions for series expansion
    %----------------------------------------------------------------------
    
    % Verify whether condition for series expansion 1 is satisfied, i.e.
    % that |varsigma-pi_{i}|<|pi_{i}-pi_{j}|, cond1=1 if satisfied.
    pi_rep = repmat(pfe.poles.',[1,length(varsigma)]);
    varsigma_rep = repmat(varsigma,[length(pfe.poles),1]);
    varsigma_pi_diff = abs(varsigma_rep-pi_rep);
    [c1_1, min_idx] = min(varsigma_pi_diff(:));
    [idx_p,idx_vs] = ind2sub(size(varsigma_pi_diff),min_idx);
    clear pi_rep; clear varsigma_rep; clear varsigma_pi_diff;
    aux = pfe.poles;
    aux(idx_p) = [];
    c1_2 = min(abs(aux-pfe.poles(idx_p)));
    clear aux;
    % Condition flag
    cond1 = (c1_1<c1_2);
    if (cond1~=1)
        stop1
    end
    
    % Verify whether condition for series expansion 2 is satisfied, i.e.
    % that |varsigma|<|pi_{i}| , cond2=1 if satisfied.
    c2_1 = abs(varsigma(end));
    c2_2 = min(abs(pfe.poles));
    % Condition flag
    cond2 = (c2_1<c2_2);
    if (cond2~=1)
        stop2
    end
    
    %----------------------------------------------------------------------
    % Series expansion approximation
    %----------------------------------------------------------------------
    approx_type = 1;
    
    if approx_type==1    
        % Highest order in series expansion, i.e. remainder of n_max+1
        n_max = R;  
        % Approximation with series expansion 1
        G_loc = rat_approx1(varsigma,pfe,n_max); 
    end
    G_loc1 = G_loc;
    
    approx_type = 2;
    if approx_type==2;
        % Highest order in series expansion, i.e. remainder of n_max+1
        n_max = 2*R;
        % Approximation with series expansion 2
        [G_loc, auxilaries] = rat_approx2(varsigma,pfe,n_max);
    end
    G_loc2 = G_loc;
    
    % Collect appproximation errors and frequency resolution
    G_err1(ii) = G_loc1(1)-G0_N(fn1/fs*N+1-nn);
    G_err2(ii) = G_loc2(2)-G0_N(fn1/fs*N+1-nn);
    dsigma1(ii) = 2*nn*df;
    dsigma2(ii) = abs(c1_1);
end

%%
%--------------------------------------------------------------------------
% Plot function and appproximation at s = j*omega0
%--------------------------------------------------------------------------
freq_loc = fn1-(1i)*varsigma;
figure(1)
plot(freq_N,db(G0_N),'k')
hold on
plot(freq_loc,db(G_loc1),'r+')
hold on
plot(freq_loc,db(G_loc2),'b+')
axis([49.6 50.6 27 32])

%--------------------------------------------------------------------------
% Study appproximation errors at s = j*omega0
%--------------------------------------------------------------------------

% Plot error G0(varsigma)-G(varsigma) at varsigma=j*nn*2*pi*df in function
% of df on a log-log scale.
idx_l = 1;
idx_r = ii;

figure(2)
plot(log(dsigma1(1:idx_r)),log(abs(G_err1(1:idx_r))),'r.-')
hold on
plot(log(dsigma1(1:idx_r)),log(abs(G_err2(1:idx_r))),'b.-')

% Plot error error G0(varsigma)-G(varsigma) at varsigma=j*nn*2*pi*df in
% function of (varsigma-pi_{min}) on a log-log scale.
hold on
plot(log(dsigma2(1:idx_r)),log(abs(G_err1(1:idx_r))),'r.-')
hold on
plot(log(dsigma2(1:idx_r)),log(abs(G_err2(1:idx_r))),'b.-')

% Some slope calculations
rico1 = log(abs(G_err1(idx_l)))-log(abs(G_err1(idx_r)))/(log(dsigma1(idx_l))-log(dsigma1(idx_r)))
rico2 = log(abs(G_err2(idx_l)))-log(abs(G_err2(idx_r)))/(log(dsigma1(idx_l))-log(dsigma1(idx_r)))

%--------------------------------------------------------------------------
% Study appproximation errors at s = j*omega0
%--------------------------------------------------------------------------







