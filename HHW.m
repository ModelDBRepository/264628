% Code by Ehsan Mirzakhalili 
% mirzakh@umich.edu, ehsan.mirzakhalili@gmail.com
% https://doi.org/10.3389/fncir.2017.00038

%% Constructing the function for passing to MATLAB ode45
% t is time
% X are the states
% W is the netwrok connectivity
% Wij=1 means that the presynaptic neuron j and the postsynaptic neuron i
% are connected
% Iext is the external intial drive that is turend off after 100 ms
%%
function F=HHW(t,X,W,Iext)
%%
N=length(W);% Number of neurons
%% Neurons' paramters (doi: 10.1371/journal.pcbi.1002062)
ENa=55;% Sodium reversal potential
EK=-90;% Potassium reversal potential
EL=-60;% Leak reversal potential
gNa=24;% Sodium conductance
gK=3;% Potassium conductance
gL=0.02;% Leak conductance
CM=1;% Capacitance
Idrive=-0.137; % Current to keep cells quiescence
%% Syanpses' parameters (Destexhe et al., 1998)
EA=0; % AMPA reversal potential
gA=1/N; % Conductance is scaled by the newtwork size. Does not apply to this study
aA=1.1; % forward rate
bA=0.19;% backward rate
TMax=1;% Maximum concentration of released neurotransmitters
VT=2;Kp=5; % constantsfor the steepness and half-activation of neurotransmitter release
%% States
%Initializing
v=X(1:N); % Membrane voltage
h=X(1*N+1:2*N);% Sodium inactivation
n=X(2*N+1:3*N);% potassium activation
s=X(3*N+1:4*N);% fraction of open receptors
%%
%Turning off the external current after 100 ms
if t>100
    Iext=0;
end
%% Steady-state and time constants
Minf=1./(1+exp((-v-30)/9.5)); % Sodium activation steady-state

Hinf=1./(1+exp((v+53)/7)); % Sodium inactivation steady-state
Th=0.37+2.78./(1+exp((v+40.5)/6));% Sodium activation time constant

Ninf=1./(1+exp((-v-30)/10));% Potassium activation steady-state
Tn=0.37+1.85./(1+exp((v+27)/15));% Potassium activation time constant
%%Currents
INa=gNa.*Minf.^3.*h.*(v-ENa);% Sodium current
IK=gK.*n.^4.*(v-EK);% Potassium current
IL=gL.*(v-EL);%Leak current
%Synapse
T=TMax./(1+exp(-(v-VT)./Kp));% neurotransmitters release
ISyn=-gA.*(W*s).*(v-EA);% AMPA current
%% Derivatives
F(1:N      ,1)=1./CM*(-INa-IK-IL+Idrive+Iext+ISyn);% Current balance
F(1*N+1:2*N,1)=(Hinf-h)./Th;% h-ODE
F(2*N+1:3*N,1)=(Ninf-n)./Tn;% n-ODE
F(3*N+1:4*N,1)=aA.*T.*(1-s)-bA.*s;%s-ODE
end