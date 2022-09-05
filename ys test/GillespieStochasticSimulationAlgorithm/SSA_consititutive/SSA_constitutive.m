%%%Daniel Charlebois - May 2011 - MATLAB v7.11 (R2010b)
%Stochastic Simulation Algorithm (SSA) - Gillespie's Direct Method (Gillespie,J.Phys.Chem.,1977).
%Default reactions and analytical solutions correspond to a simple model of
%gene expression for a single cell (Kaern et al.,Nat.Rev.Genet.,2005): 
%Pro-->M-->P, with mRNA (M-->0) and protein (P-->0) decay, when expression 
%is assumed to be constitutive and cell volume is constant and set to unity.

clc; clear all; close all;

%% initialization
%%%general simulation parameters
t=0; %start time 
t_end=1000; %end time
t_sample=1; %sample interval for gathering data
k=1; %counter for waitbar update 
alpha=10^2; %parameter for updating waitbar (increase for shorter runtime)

%%%model parameters
%initial numbers for each chemical species 
Pro = 1; %promoter
M = 5; %mRNA
P = 100; %protein
%rate constants 
kM = 5.5; %mRNA production (transcription) 
kP = 1; %protein production (translation)
dM = 0.5; %mRNA decay
dP = 0.05; %protein decay

%%%arrays to store results
j=1; %counter for arrays
t_array(1,t_end/t_sample+1)=0; t_array(1,j)=t; %time array and initial value
Pro_array(1,t_end/t_sample+1)=0; Pro_array(1,j)=Pro; %promoter array and initial value
M_array(1,t_end/t_sample+1)=0; M_array(1,j)=M; %mRNA array and initial value
P_array(1,t_end/t_sample+1)=0; P_array(1,j)=P; %protein array and initial value

%% SSA
tic %start timing the Gillespie loop
w=waitbar(0,'running SSA...');
while t < t_end,
    
    %calculate rxn propensities
    h = [kM*Pro kP*M dM*M dP*P];
    %combined rxn hazard
    h0 = sum(h);
    
    %calculate time to next event
    r1=rand;
    while r1 == 0,
        r1=rand;
    end
    t_next = ((1/h0)*(log(1/r1)));
    
    %update time
    t = t + t_next;
    
    %determine next reaction
    i=1; mu=0; amu=0; r2=rand;
    while amu < r2*h0,
        mu = mu + 1;
        amu = amu + h(i); 
        i = i + 1;
    end
    
    %reactions
    if mu == 1 %transcription
        M = M + 1;
    elseif mu == 2 %translation
        P = P + 1;
    elseif mu == 3 %mRNA decay
        M = M - 1;
    elseif mu == 4 %protein decay
        P = P - 1;
    end
    
    %store/output time and species
    if t >= j*t_sample
        j=j+1;
        t_array(1,j)=j;
        Pro_array(1,j)=Pro;
        M_array(1,j)=M;
        P_array(1,j)=P;
    end    
    
    %update waitbar
    if t >= k*alpha*t_sample
        k=k+1;
        waitbar(t/t_end)
    end
end 
close(w)
toc

%% analytical results for model
mean_mRNA_theoretical = kM/dM %steady-state mean mRNA
noise_mRNA_theoretical = sqrt(1/mean_mRNA_theoretical) %steady-state noise mRNA
mean_protein_theoretical = (kM*kP)/(dM*dP) %steady-state mean protein
noise_protein_theoretical = sqrt((1/mean_protein_theoretical)+ (1/(mean_mRNA_theoretical*(1+(dM/dP))))) %steady-state noise protein

%% output simulation statistics
mean_mRNA_simulation = mean(M_array)
stdev_mRNA_simulation = std(M_array);
noise_mRNA_simulation = stdev_mRNA_simulation/mean_mRNA_simulation
mean_protein_simulation = mean(P_array) 
stdev_protein_simulation = std(P_array);
noise_protein_simulation = stdev_protein_simulation/mean_protein_simulation

%% plots
figure; 
subplot(2,2,1);
plot(t_array,M_array); %mRNA time series
xlabel('time (s)');
ylabel('mRNA no.');
subplot(2,2,2);
hist(M_array); %mRNA histogram
xlabel('mRNA no.');
ylabel('counts');
subplot(2,2,3);
plot(t_array,P_array); %protein time series
xlabel('time (s)');
ylabel('protein no.');
subplot(2,2,4);
hist(P_array); %protein histogram
xlabel('protein no.');
ylabel('counts');