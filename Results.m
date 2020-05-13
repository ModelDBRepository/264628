% Code by Ehsan Mirzakhalili 
% mirzakh@umich.edu, ehsan.mirzakhalili@gmail.com
% https://doi.org/10.3389/fncir.2017.00038
clc;
clear;
rng(1); % Fixing the random number seed for reproducibility
TEnd=4000; % Lenght of the simulation
Vth=0; % Threshold for action potential detection
V0=-65.11; % Resting potential
N=200; % Number of neurons
IRand=rand(N,1);% Intial random drive
Hinf=1./(1+exp((V0+53)/7)); % Steady-state of sodium inactivation at rest
Ninf=1./(1+exp((-V0-30)/10));% Steady-state of potassium activation at rest
X0=[V0*ones(1,N) Hinf*ones(1,N) Ninf*ones(1,N) 0*ones(1,N)];%initializing the states
opts = odeset('Abstol',1e-2);% ODE options
Mean=[10 30];% Modes of the degree distribution
Weights=[0.5 0.5];% Weight of each mode
WO=MultiDistribution(N,Mean,Weights); % Constructing the unimpaired/original network;
NonZero=find(WO);%Finding the nonzero elements of the network
L=length(NonZero);% number of synapses (nonzero elements)
[TO,XO] = ode45(@(t,X) HHW(t,X,WO,IRand),[0 TEnd],X0,opts);% Calling ODE45 to solve the equations
VO=XO(:,1:N);% Taking the voltages from the solved states
SpikesO=Raster(VO,TO,Vth);% Detecting spikes
%% Random Impairemnt (Scenario 1)
jj=0.8;% Percentage of impairment
k=0.25;% Level of impairment
r=randperm(L,floor(L*jj));% Randomly selecting synapses to be impaired
W1=WO;
W1(NonZero(r))=1-k; % Performing the impairment randomly
[T1,X1] = ode45(@(t,X) HHW(t,X,W1,IRand),[0 TEnd],X0,opts);
V1=X1(:,1:N);
Spikes1=Raster(V1,T1,Vth);
%% Impairment based on outdegree (Scenario 2)
% Finding outdegrees. Note for Wij=1  
% the presynaptic neuron is j and the postsynaptic neuron is i
Dout=sum(WO,1);
L=sum(Dout);
% Sorting degrees from highest to lowest
[~,iD]=sort(Dout,'descend');
% Finding synapses with highest outdegree
Synapses=[];
for j=1:N
    Synapses=[Synapses;find(WO(:,iD(j)))+(iD(j)-1)*N];
end
jj=0.8;% Percentage of impairment
k=0.25;% Level of impairment
W2=WO;
W2(Synapses(1:floor(L*jj)))=1-k;% Performing the impairment based on the outdegree
[T2,X2] = ode45(@(t,X) HHW(t,X,W2,IRand),[0 TEnd],X0,opts);
V2=X2(:,1:N);
Spikes2=Raster(V2,T2,Vth);
%% Impairment based on the activity (Scenario 3)
% Calculating activity of each neuorn for the unimpaired network
for j=1:N
    SpikeCounts(j)=sum(SpikesO(:,2)==j);
end
% Sorting activity from highest to lowest
[~,iD]=sort(SpikeCounts,'descend');
Synapses=[];
% Finding synapses of neurons with highest activity
for j=1:N
    Synapses=[Synapses;find(WO(:,iD(j)))+(iD(j)-1)*N];
end
jj=0.8;% Percentage of impairment
k=0.25;% Level of impairment
W3=WO;
W3(Synapses(1:floor(L*jj)))=1-k;% Performing the impairment based on activity
[T3,X3] = ode45(@(t,X) HHW(t,X,W3,IRand),[0 TEnd],X0,opts);
V3=X3(:,1:N);
Spikes3=Raster(V3,T3,Vth);
%%
subplot(2,2,1);plot(SpikesO(:,1),SpikesO(:,2),'.');
xlabel('Time [ms]');ylabel('Neuron #');
xlim([0,4000]);
title('Persistent activity: no impairment')
%
subplot(2,2,2);plot(Spikes1(:,1),Spikes1(:,2),'.');
xlabel('Time [ms]');ylabel('Neuron #');
title('Scenario 1: random impairment')
xlim([0,4000]);
%
subplot(2,2,3);plot(Spikes2(:,1),Spikes2(:,2),'.');
xlabel('Time [ms]');ylabel('Neuron #');
title('Scenario 2: outdegree based impairment')
xlim([0,4000]);
%
subplot(2,2,4);plot(Spikes3(:,1),Spikes3(:,2),'.');
xlabel('Time [ms]');ylabel('Neuron #');
title('Scenario 3: Activity based impairment')
xlim([0,4000]);