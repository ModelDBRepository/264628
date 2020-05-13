% Code by Ehsan Mirzakhalili 
% mirzakh@umich.edu, ehsan.mirzakhalili@gmail.com
% https://doi.org/10.3389/fncir.2017.00038

%% This functions constructs the network
% W is the connectivity matrix
% N is the number of nodes (neurons)
% Mean are the modes of the distributions
% Weights are the weights for each mode in the distribution

function W=MultiDistribution(N,Mean,Weights)
Num=length(Mean); %% Finding number of modes
Min=0;% The network needs to be connected
Max=N;% The largest degree connot be larger than the network size
Width=1;% The resolution is 1 degree (the netwrok is binary)
Bins=(Max-Min)/Width;
D=[];A=0;
%% Creading the total degree distribution. The distribution is a combination of independent Poission distributions
for i=1:Bins
    x1=Min+(i-1)*Width;
    x2=Min+(i  )*Width;
    y1=0;y2=0;
    for j=1:Num
        pd=makedist('pois',Mean(j));
        y1=y1+pdf(pd,x1)*Weights(j);
        y2=y2+pdf(pd,x2)*Weights(j);
    end
    A=A+Width*(y1+y2)/2;
    C=Width*(y1+y2)/2*N;
    C=round(C);
    D=[D,randi([x1,x2],[1,C])];
end
%% Making sure the number of nodes match the desired network size
while length(D)<N
    D=[D,D(randperm(length(D),1))];
end
D=D(randperm(length(D),N));
%% Creating List of nodes randomly
% Each node is repeated by its degree
R=randperm(N,N);
Degree=zeros(1,N);
List=[];
for i=1:N
    Degree(i)=D(R(i));
    List=[List,i*ones(1,Degree(i))];
end
NL=length(List);
List=List(randperm(NL,NL));
if mod(NL,2)==1
    List(randperm(NL,1))=[];
    NL=NL-1;
end
W=zeros(N);
%% Creating Connectivity
Failed=0;
while isempty(List)==0
    %% Connecting the nodes randomly and removing them from the list
    % the networks with self-connections are rejected
    R1=randperm(NL,1);
    R2=randperm(NL,1);
    if W(List(R1),List(R2))~=1 && List(R1)~=List(R2)
        W(List(R1),List(R2))=1;
        if R1>R2
            List(R1)=[];
            List(R2)=[];
        else
            List(R2)=[];
            List(R1)=[];
        end
        NL=NL-2;
        Failed=0;
    else
        Failed=Failed+1;
        if Failed>1000
            break;
        end
    end
end