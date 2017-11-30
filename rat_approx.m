clear all,close all,clc

% dsp parameters
N = 4096;
fs = 256;
df = fs/N;
freq = [1:N/2-1]*df;
Omega = (1i)*2*pi*freq;

% System 
%%%%%%%%%

% first second-order system 
G_dc1 = 2;
fn1 = 50;
varsigma1 = 0.005;
omega_n1 = 2*pi*fn1;

b1 = -(1i)*G_dc1*omega_n1/(2*sqrt(1-varsigma1^2));
p1 = -varsigma1*omega_n1+(1i)*omega_n1*sqrt(1-varsigma1^2);
G1 = b1./(Omega-p1) + conj(b1)./(Omega-conj(p1));

% another second-order system
G_dc2 = 2;
fn2 = 45;
varsigma2 = 0.05;
omega_n2 = 2*pi*fn2;

b2 = -(1i)*G_dc2*omega_n2/(2*sqrt(1-varsigma2^2));
p2 = -varsigma2*omega_n2+(1i)*omega_n2*sqrt(1-varsigma2^2);
G2 = b2./(Omega-p2) + conj(b2)./(Omega-conj(p2));

G = G1+G2;

% Plot system response
figure(1)
semilogx(freq,db(G),freq,db(G1))

% Collect poles and constant coefficients
poles = [p1,conj(p1);p2,conj(p2)];
coef = [b1,conj(b1);b2,conj(b2)];


%% Linearized least squares approximation
R = 1;
fc = 50;    % frequency of interest and center frequency of local band




