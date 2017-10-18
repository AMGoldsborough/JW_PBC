function JW_PBC(L,Jstr,Jdis,Pdist,Jseed)
%JW_PBC(L,Jstr,Jdis,Pdist,Jseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JW_PBC
% ED by Jordan Wigner for XX model with PBC
% based on 10.1103/PhysRevB.57.11457
%
% Andrew Goldsborough - 21/02/2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%inputs
% L = 10;         %chain length
% Jstr = 1;       %overall J strength
% Jdis = 2;       %disorder strength
% Pdist = 4;      %coupling distribution
% Jseed = 2;      %seed for rng, 0 => shuffle

%coupling distribution
%0 => manual
%1 => 2 theta(K-1/2)
%2 => 1
%3 => uniform around Jstr normalised by Jstr
%4 => uniform around Jstr un-normalised Jdis
%5 => box distribution of Hikihara AF (10.1103/PhysRevB.60.12116)
%6 => Laflorencie's infinite disorder distribution (10.1103/PhysRevB.72.140408)

if mod(L,2) == 1
    error('odd lattice sizes not currently implemented');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate couplings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set seed, set by clock if zero
if Jseed > 0
    rng(Jseed);
else
    rng('shuffle');
end

%set the probability distribution
if Pdist==0
    %custom set
    J = [0.1,Jdis,0.1,1,Jdis,1,0.1,Jdis,0.1,0.01];
elseif Pdist==1
    %P(K) = 2 theta(K-1/2)
    J = zeros(1,L) + Jstr*random('unif',0.5,1,[1,L]);
elseif Pdist==2
    %P(K) = 1
    J = zeros(1,L) + Jstr*(rand(L,1));
elseif Pdist==3
    %uniform around Jstr normalised by Jstr
    J = zeros(1,L) + Jstr + Jstr*Jdis*(rand(1,L) - 0.5);
elseif Pdist==4
    %uniform around Jstr un-normalised Jdis
    J = zeros(1,L) + Jstr + Jdis*(rand(1,L) - 0.5);
elseif Pdist==5
    %box distribution of Hikihara AF
    J = zeros(1,L) + Jdis*random('unif',0,1,[1,L]);
elseif Pdist==6
    %Laflorencie's infinite disorder distribution
    J = rand(1,L).^Jdis;
end

rng('shuffle');

%print inputs
fprintf('L = %d, Jstr = %f, Jdis = %f, Pdist = %d, Jseed = %d\n',L,Jstr,Jdis,Pdist,Jseed);

%print J to file
fprintf('printing interaction strengths\n');
fname = strcat('./J/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_J.txt');
fidJ = fopen(fname, 'w');
for i = 1:L
    fprintf(fidJ,'%.15e\n',J(i));
end
fclose(fidJ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagonalise by JW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make Hamiltonian matrix
H = zeros(L,L);
H(1,2) = J(1)/2;
for i = 2:L-1
    H(i,i-1) = J(i-1)/2;
    H(i,i+1) = J(i)/2;
end
H(L,L-1) = J(L-1)/2;

%different BCs when number of particels is odd or even (note half filling)
if mod(L/2,2) == 1
    %PBC
    H(1,L) = J(L)/2;
    H(L,1) = J(L)/2;
else %mod(L/2,2) == 0
    %aPBC
    H(1,L) = -J(L)/2;
    H(L,1) = -J(L)/2;
end

%diagonalise
[v,nrg] = eig(H,'vector');

%ground state is found by summing the eigs up to the filling fraction, L/2
energy = sum(nrg(1:L/2));

%print energy to file
fprintf('printing energy\n');
fname = strcat('./energy/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_energy_JW_PBC.txt');
fidenergy = fopen(fname, 'w');
fprintf(fidenergy,'%.15e\n',energy);
fclose(fidenergy);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correlation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%matrix of <Ci'.Cj>
Cm = v(:,1:L/2) * v(:,1:L/2).';

%matrix of <Ai.Bj>
ABm = eye(L) - 2*Cm;

%Sz corr
fprintf('printing Sz correlation functions\n');
fname = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_Szcorr_JW_PBC.txt');
fidSzcorr = fopen(fname, 'w');

%Sx corr
fprintf('printing Sx correlation functions\n');
fname = strcat('./Sxcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_Sxcorr_JW_PBC.txt');
fidSxcorr = fopen(fname, 'w');

for i = 1:(L-1)
    for j = (i+1):L
        Sz_corr = 0.25*det([ABm(i,i) ABm(i,j); ABm(j,i) ABm(j,j)]);
        fprintf(fidSzcorr,'%d %d %.15e\n',i,j,Sz_corr);
        
        Sx_corr = 0.25*det(-ABm(i:(j-1),(i+1):j));
        fprintf(fidSxcorr,'%d %d %.15e\n',i,j,Sx_corr);
    end
end

fclose(fidSzcorr);
fclose(fidSxcorr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%entanglement entropy (see 10.1103/PhysRevB.72.140408)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('printing ee\n');
fname = strcat('./ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_ee_JW_PBC.txt');
fidee = fopen(fname, 'w');

%set normalisation
ee_eigs = 0.5;
norm = -sum(ee_eigs.*log(ee_eigs) + (1-ee_eigs).*log(1-ee_eigs));

for i = 2:L
    for j = i:L
        ee_eigs = eig(0.5*(Cm(i:j,i:j) + Cm(i:j,i:j)'));
        
        %sort stability
        ee_eigs(ee_eigs>=1) = 1-1e-16;
        ee_eigs(ee_eigs<=0) = 1e-16;

        %calculate ee and normalise
        ee = -sum(ee_eigs.*log(ee_eigs) + (1-ee_eigs).*log(1-ee_eigs));
        fprintf(fidee,'%d %d %.15e\n',i-1,j,ee/norm);
    end
end

fclose(fidee);
toc
end



