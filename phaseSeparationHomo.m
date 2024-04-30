%% homo interactions -- total energy calculation
% A + A <--> AA
%constants

kd = 1e-2; % dimer binding affinity
A_tot = 0.1e-3; % total concentration
tot_protein = A_tot;

% values in the absence of condensates
AA_original = (1/8)*(4*A_tot+kd-sqrt((-4*A_tot-kd)^2-16*A_tot^2));
A_original = A_tot - 2*AA_original;

%%%
c_sat = 0.1; % saturation concentration - limit of concentration in dense phase
chemPot_constant = 1e-1; %constant to induce energetic cost a dense phase concentrations
                         %approaching saturation concentration

%variables
part = [1:1:200]; %range of partition coefficients
Vf_log = [-4:0.1:0]; %range of volume fractions - log10
Vf = 10.^Vf_log;

[X,Y] = meshgrid(Vf_log,part); %xy for plotting

% preallocate
A_tot_cond = zeros(length(part),length(Vf_log));
A_free_cond = zeros(length(part),length(Vf_log));

AA_cond = zeros(length(part),length(Vf_log));

A_tot_dil = zeros(length(part),length(Vf_log));
A_free_dil = zeros(length(part),length(Vf_log));

AA_dil = zeros(length(part),length(Vf_log));

A_free_tot = zeros(length(part),length(Vf_log));
AA_tot = zeros(length(part),length(Vf_log));

S_bound = zeros(length(part),length(Vf_log));
S_A_free = zeros(length(part),length(Vf_log));

cond_tot_protein = zeros(length(part),length(Vf_log));
u_cond = zeros(length(part),length(Vf_log));

%calulations at varying partition coefficient and volume fraction
for i=1:length(part)
    for j=1:length(Vf_log)

        % calculate condensate concentrations
        A_tot_cond(i,j) = A_tot/(Vf(j)+(1-Vf(j))/part(i));

        AA_cond(i,j) = (1/8)*(4*A_tot_cond(i,j)+kd-sqrt((-4*A_tot_cond(i,j)-kd)^2-16*A_tot_cond(i,j)^2));

        A_free_cond(i,j) = A_tot_cond(i,j)-2*AA_cond(i,j);

        %calculate dilute region concentrations
        A_tot_dil(i,j) = A_tot/(Vf(j)*(part(i)-1)+1);

        AA_dil(i,j) = (1/8)*(4*A_tot_dil(i,j)+kd-sqrt((-4*A_tot_dil(i,j)-kd)^2-16*A_tot_dil(i,j)^2));

        A_free_dil(i,j) = A_tot_dil(i,j)-2*AA_dil(i,j);

        cond_tot_protein(i,j) = A_tot_cond(i,j);

        %calculate chemical potential at high concentration dense phase
        if (c_sat/(c_sat-cond_tot_protein(i,j)))>0
            u_cond(i,j) = chemPot_constant*Vf(j)*(exp((c_sat/(c_sat-cond_tot_protein(i,j)))^2 -1) - exp((c_sat/(c_sat-tot_protein))^2 -1)); %
        else
            u_cond(i,j)=nan;
        end

        %calculate total concentrations
        A_free_tot(i,j) = Vf(j)*A_free_cond(i,j)+(1-Vf(j))*A_free_dil(i,j);
        AA_tot(i,j) = Vf(j)*AA_cond(i,j)+(1-Vf(j))*AA_dil(i,j);

        %calculate entropy
        S_bound(i,j) = AA_cond(i,j)*Vf(j)*log(AA_cond(i,j))+AA_dil(i,j)*(1-Vf(j))*log(AA_dil(i,j))-AA_original*log(AA_original);
        S_A_free(i,j) = A_free_cond(i,j)*Vf(j)*log(A_free_cond(i,j))+A_free_dil(i,j)*(1-Vf(j))*log(A_free_dil(i,j))-A_original*log(A_original);

    end
end

%calculate free energy of binding
kd_eff = A_free_tot.^2./AA_tot;
dG_bind_permole = log(kd_eff/kd);
dG_bind_tot = dG_bind_permole*tot_protein;

%sum entropy
S_tot = S_A_free+S_bound;
S_tot_permole = S_tot/tot_protein;

%calculate total condensate concentration
cond_tot_protein = A_tot_cond;
dil_tot_protein = A_tot_dil;

u_cond(u_cond>0.0002) = nan; % remove large energies for plotting clarity

%total energy
E_tot = dG_bind_tot+S_tot + u_cond;
E_tot_permole = dG_bind_permole+S_tot_permole;

conc_kd = tot_protein/kd;
minEnergy = min(E_tot(:));
[p,v] = find(E_tot == minEnergy);
part_min = part(p);
vf_min = Vf(v);
condConc = cond_tot_protein(p,v);


figure()
surf(X,Y,E_tot)
figure();contourf(X,Y,S_tot,20)
