%% hetero interactions -- total energy calculation
% A + B <--> AB
%constants

kd = 1e-3; % dimer binding affinity
A_tot = 0.05e-4; % total concentration of A
B_tot = A_tot; % total concentration of B
tot_protein = A_tot + B_tot; % total combined concentration

B_mult = 1; %number of A that can bind one B molecule

% values in the absence of condensates
AB_original = 0.5*(A_tot+B_tot+kd-sqrt((-A_tot-B_tot-kd)^2-4*A_tot*B_tot));
A_original = A_tot - AB_original;
B_original = B_tot - AB_original;

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

B_tot_cond = zeros(length(part),length(Vf_log));
B_free_cond = zeros(length(part),length(Vf_log));

AB_cond = zeros(length(part),length(Vf_log));

A_tot_dil = zeros(length(part),length(Vf_log));
A_free_dil = zeros(length(part),length(Vf_log));

B_tot_dil = zeros(length(part),length(Vf_log));
B_free_dil = zeros(length(part),length(Vf_log));

AB_dil = zeros(length(part),length(Vf_log));

A_free_tot = zeros(length(part),length(Vf_log));
B_free_tot = zeros(length(part),length(Vf_log));
AB_tot = zeros(length(part),length(Vf_log));

S_bound = zeros(length(part),length(Vf_log));
S_A_free = zeros(length(part),length(Vf_log));
S_B_free = zeros(length(part),length(Vf_log));

cond_tot_protein = zeros(length(part),length(Vf_log));
u_cond = zeros(length(part),length(Vf_log));

%%%
%calulations at varying partition coefficient and volume fraction
for i=1:length(part)
    for j=1:length(Vf_log)

        % calculate condensate concentrations
        A_tot_cond(i,j) = A_tot/(Vf(j)+(1-Vf(j))/part(i));
        B_tot_cond(i,j) = B_tot/(Vf(j)+(1-Vf(j))/part(i));

        AB_cond(i,j) = 0.5*(A_tot_cond(i,j)+B_tot_cond(i,j)+kd-sqrt((-A_tot_cond(i,j)-B_tot_cond(i,j)-kd)^2-4*A_tot_cond(i,j)*B_tot_cond(i,j)));

        A_free_cond(i,j) = A_tot_cond(i,j)-AB_cond(i,j);
        B_free_cond(i,j) = B_tot_cond(i,j)-AB_cond(i,j);

        %calculate dilute region concentrations
        A_tot_dil(i,j) = A_tot/(Vf(j)*(part(i)-1)+1);
        B_tot_dil(i,j) = B_tot/(Vf(j)*(part(i)-1)+1);

        AB_dil(i,j) = 0.5*(A_tot_dil(i,j)+B_tot_dil(i,j)+kd-sqrt((-A_tot_dil(i,j)-B_tot_dil(i,j)-kd)^2-4*A_tot_dil(i,j)*B_tot_dil(i,j)));

        A_free_dil(i,j) = A_tot_dil(i,j)-AB_dil(i,j);
        B_free_dil(i,j) = B_tot_dil(i,j)-AB_dil(i,j);

        cond_tot_protein(i,j) = A_tot_cond(i,j)+B_tot_cond(i,j);

        %%%%
        if (c_sat/(c_sat-cond_tot_protein(i,j)))>0
            u_cond(i,j) = chemPot_constant*Vf(j)*(exp((c_sat/(c_sat-cond_tot_protein(i,j)))^2 -1) - exp((c_sat/(c_sat-tot_protein))^2 -1)); %
        else
            u_cond(i,j)=nan;
        end

        %calculate total concentrations
        A_free_tot(i,j) = Vf(j)*A_free_cond(i,j)+(1-Vf(j))*A_free_dil(i,j);
        B_free_tot(i,j) = Vf(j)*B_free_cond(i,j)+(1-Vf(j))*B_free_dil(i,j);
        AB_tot(i,j) = Vf(j)*AB_cond(i,j)+(1-Vf(j))*AB_dil(i,j);

        %calculate entropy
        S_bound(i,j) = AB_cond(i,j)*Vf(j)*log(AB_cond(i,j))+AB_dil(i,j)*(1-Vf(j))*log(AB_dil(i,j))-AB_original*log(AB_original);
        S_A_free(i,j) = A_free_cond(i,j)*Vf(j)*log(A_free_cond(i,j))+A_free_dil(i,j)*(1-Vf(j))*log(A_free_dil(i,j))-A_original*log(A_original);
        S_B_free(i,j) = B_free_cond(i,j)/B_mult*Vf(j)*log(B_free_cond(i,j)/B_mult)+B_free_dil(i,j)/B_mult*(1-Vf(j))*log(B_free_dil(i,j)/B_mult)-B_original/B_mult*log(B_original/B_mult);

    end
end

%calculate free energy of binding
kd_eff = A_free_tot.*B_free_tot./AB_tot;
dG_bind_permole = log(kd_eff/kd);
dG_bind_tot = dG_bind_permole*tot_protein;

%sum entropy
S_tot = S_A_free+S_B_free+S_bound;
S_tot_permole = S_tot/tot_protein;

%calculate total condensate concentration
cond_tot_protein = A_tot_cond+B_tot_cond;
dil_tot_protein = A_tot_dil+B_tot_dil;

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
figure();contourf(X,Y,E_tot,20)


