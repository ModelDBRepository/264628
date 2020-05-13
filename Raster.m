% Code by Ehsan Mirzakhalili 
% mirzakh@umich.edu, ehsan.mirzakhalili@gmail.com
% https://doi.org/10.3389/fncir.2017.00038
%% Finding peaks to detect action potentials
% For voltages that are larger than Vth, we check if V of neuron j at time i is a local
% maximum
function Spikes=Raster(V,T,Vth)
N=size(V,2);
S=size(V,1);
J=0;
for i=2:S-1
    for j=1:N
        if V(i,j)>Vth
            if (V(i,j)-V(i-1,j))*(V(i,j)-V(i+1,j))>0
                J=J+1;
                Spikes(J,1)=T(i);
                Spikes(J,2)=j;
            end
        end
    end
end