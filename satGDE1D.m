function [UI, SOL] = satGDE1D

% CONSTANTS
M_w = 18e-3; % [kg/mol] molar mass of water
rho_w = 997; % [kg/m^3] density of water at reference temperature
T_0 = 273.15; % [K] zero degrees Celsius
R = 8.31446; % [J*mol/K] universal gas constant
F = 96485.333; % [C/mol] Faraday constant
P_ref = 101325; % [Pa] reference pressure
T_ref = T_0+25; % [K] reference temperature

% OPERATING CONDITIONS
P = 1e5; % [Pa] total pressure in gas channel
RH = 0.9; % [-] relative humidity in gas channel
p_C_GC = 0; % [Pa] capillary pressure at GDL/GC interface
p_C_LC = 0; % [Pa] capillary pressure at Liquid channel interface
T_C = T_0+30; % [K] temperature of gas channel
T_L = T_0+30; % [K] temperature of liquid channel
alpha_CO2 = 0.9999; % [-] mole fraction of carbon dioxide in dry feed gas
alpha_CO = 0.00005; % [-] mole fraction of carbon monoxide in dry feed gas
alpha_H2 = 0.00005; % [-] mole fraction of hydrogen in dry feed gas
C_I = 100; % [mol/m^3] ionic strength of electrolyte
U = [0:-0.025:-0.75]; % [V] applied voltage at current collector

% ELECTROCHEMICAL PARAMETERS
i_0_COER = @(T) 7.25e9.*exp(-100e3/R./T); % [A/m^2] exchange current density of COER (alkaline)
i_0_HER = @(T,pH) 8.84e7.*exp(-(83+1.*pH)*1e3/R./T); % [A/m^2] exchange current density of HER (alkaline)
beta_COER = 1; % [-] COER symmetry factor
beta_HER = 0.44; % [-] HER symmetry factor
U_0_COER = -0.11; % [V] COER standard electrode potential relative to SHE
U_0_HER = 0; % [V] HER standard electrode potential relative to SHE
C_ref_CO2 = 1e3; % [mol/m^3] reference aqeuous concentraiton for CO2
p_ref_H2 = 1; % [atm] reference pressure for H2
gamma_BV_H = 1; % Proton correction exponent for concentration dependent Butler-Volmer eq
gamma_BV_CO2 = 1.5; % CO2 correction exponent for concentration dependent Butler-Volmer eq
Pi_COER = @(T) 0.240*T/298; % [V] Peltier coefficient (irreversible loss) for COER 
Pi_HER = @(T) 0.013*T/298; % [V] Peltier coefficient (irreversible loss) for HER 

% MATERIAL PARAMETERS
L = [300 45 5]*1e-6; % [m] gas diffusion electrode domain thicccnesses
a_CCL = 1e7; % [m^2/m^3] ECSA density of CCL
H_ec = 42e3; % [J/mol] molar enthalphy of evaporation/condensation
k_GDL = 1.2; % [W/(m*K)] thermal conductivity of GDL ????????????????????????
k_MPL = 0.2; % [W/(m*K)] thermal conductivity of MPL ????????????????????????
k_CL = 0.27; % [W/(m*K)] thermal conductivity of CL ?????????????????????????
sigma_e_GDL = 1250; % [S/m] electrical conductivity of GDL
sigma_e_MPL = 1250/6; % [S/m] electrical conductivity of GDL (ASSUMED FOR NOW)
sigma_e_CL = 350; % [S/m] electrical conductivity of CL ?????????????????????
s_im_GDL = 1e-9; % [-] immobile liquid water saturation of GDL
s_im_MPL = 1e-9; % [-] immobile liquid water saturation of MPL
s_im_CL = 1e-9; % [-] immobile liquid water saturation of CL
eps_p_GDL = 0.7; % [-] porosity of GDL
eps_p_MPL = 0.3; % [-] porosity of MPL
eps_p_CL = 0.42; % [-] porosity of CL
kappa_L_GDL = 0.8e-11; % [m^2] absolute permeability of GDL
kappa_L_MPL = 5e-14; % [m^2] absolute permeability of MPL
kappa_L_CL = 1e-13; % [m^2] absolute permeability of CL
tau_GDL = 1.6; % [-] pore tortuosity of GDL
tau_MPL = 1.6; % [-] pore tortuosity of MPL
tau_CL = 1.6; % [-] pore tortuosity of CL
theta_GDL = 93; % [°] intrinsic mean contact angle of GDL
theta_MPL = 110; % [°] intrinsic mean contact angle of MPL
theta_CL = 84; % [°] intrinsic mean contact angle of CL

% WATER CONSTITUTIVE RELATIONSHIPS
P_sat_o = @(T) exp(23.1963-3816.44./(T-46.13)); % [Pa] uncorrected saturation pressure of water vapor
P_sat = @(T,P_C) P_sat_o(T).*exp(P_C.*(M_w/rho_w)./(R.*T)); % [Pa] capillary pressure corrected saturation pressure of water vapor
mu = @(T) 1e-3*exp(-3.63148+542.05./(T-144.15)); % [Pa*s] dynamic viscosity of liquid water

% DISSOLVED SPECIES CONSTITUITIVE RELATIONSHIPS
H_CO2 = @(T) (1/101325)*34.*exp(2400.*(1./T-1/298)); % [mol/m^3/Pa] Henry's coefficient for carbon dioxide
K = @(deltaH_rxn,deltaS_rxn,T) exp(deltaS_rxn/R).*exp(-deltaH_rxn./(R.*T));
deltaS_w = -80.66; % [J/mol/K] reaction entropy change for water dissociation
deltaS_1 = -96.31; % [J/mol/K] reaction entropy change for carbonation step 1 (proton/acidic form)
deltaS_2 = -148.1; % [J/mol/K] reaction entropy change for carbonation step 2 (proton/acidic form)
deltaH_w = 55.84e3; % [J/mol] reaction enthalpy change for water dissociation
deltaH_1 = 7.64e3; % [J/mol] reaction enthalpy change for carbonation step 1 (proton/acidic form)
deltaH_2 = 14.85; % [J/mol] reaction enthalpy change for carbonation step 2 (proton/acidic form)
Kw = @(T) K(deltaH_w,deltaS_w,T); % [mol^2/m^6] H2O <-> H^+ + OH^-, water dissociation
K1 = @(T) K(deltaH_1,deltaS_1,T); % [mol/m^3] CO2(aq) + H2O <-> H^+ + HCO3^-, carbonation step 1 (proton/acidic form)
K2 = @(T) K(deltaH_2,deltaS_2,T); % [mol/m^3] HCO3^- <-> H^+ + CO3^2-, carbonation step 2 (proton/acidic form)
K3 = @(T) K1(T)./Kw(T); % [m^3/mol] CO2(aq) + OH^- <-> HCO3^-, carbonation step 1 (hydroxide,/basic form)
K4 = @(T) K2(T)./Kw(T); % [m^3/mol] HCO3^- + OH^- <-> H2O + CO3^2-, carbonation step 2 (hydroxide,/basic form)
k1f =         3.71e-2; % [1/s]
k1r = @(T) k1f./K1(T); % [m^3/mol/s]
k2f =           59.44; % [1/s]
k2r = @(T) k2f./K2(T); % [m^3/mol/s]
k3f =            2.23; % [m^3/mol/s]
k3r = @(T) k3f./K3(T); % [1/s]
k4f =           6.0e6; % [m^3/mol/s]
k4r = @(T) k4f./K4(T); % [1/s]

% GASEOUS DIFFUSION COEFFICIENTS (BULK)
D_ab = @(nu_p_a,nu_p_b,M_a,M_b,P,T) ...
    1e-4*(1e-3*T.^1.75.*(1/M_a+1/M_b)^0.5)./(P/101325*(nu_p_a^0.33+nu_p_b^0.33)^2); %[m^2/s] binary diffusion coefficient
nu_p_H2O = 12.7; % [-] diffusion volume of gaseous species
nu_p_CO2 = 26.9; % [-] diffusion volume of gaseous species
nu_p_CO = 18.9; % [-] diffusion volume of gaseous species
nu_p_H2 = 7.07; % [-] diffusion volume of gaseous species
M_H2O = 18.02; % [g/mol] molar mass of gaseous species
M_CO2 = 44.01; % [g/mol] molar mass of gaseous species
M_CO = 28.01; % [g/mol] molar mass of gaseous species
M_H2 = 2.02; % [g/mol] molar mass of gaseous species 
D_H2O_CO2_ref = D_ab(nu_p_H2O,nu_p_CO2,M_H2O,M_CO2,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, H2O in CO2
D_H2O_CO_ref = D_ab(nu_p_H2O,nu_p_CO,M_H2O,M_CO,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, H2O in CO2
D_H2O_H2_ref = D_ab(nu_p_H2O,nu_p_H2,M_H2O,M_H2,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, H2O in H2
D_CO_CO2_ref = D_ab(nu_p_CO,nu_p_CO2,M_CO,M_CO2,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, CO in CO2
D_H2_CO_ref = D_ab(nu_p_CO,nu_p_H2,M_CO,M_H2,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, CO in H2
D_H2_CO2_ref = D_ab(nu_p_CO2,nu_p_H2,M_CO2,M_H2,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, CO2 in H2 

% LIQUID DIFFUSION COEFFICIENTS (BULK)
D_CO2_L = @(T) 2.17e-9*exp(-2345.*(1./T-1/303)); % [m^2/s] diffusion coefficient, CO2 in water at temperature T/[K]
D_OH = @(T) 2.89e-9*exp(-1750.*(1./T-1/273)); % [m^2/s] diffusion coefficient, OH^- in water at temperature T/[K]
D_HCO3 = @(T) 7.016e-9*(T/204.03 - 1); % [m^2/s] diffusion coefficient, HCO3^- in water at temperature T/[K]
D_CO3 = @(T) 5.447e-9*(T/210.26 - 1); % [m^2/s] diffusion coefficient, CO3^2- in water at temperature T/[K]

% MODEL PARAMETERIZATION
BV = @(i_0,c,c_ref,gamma,beta,T,eta) i_0.*(c./c_ref).^gamma.*exp(-beta*F./(R.*T).*eta); % [A/m^2] Butler-Volmer eq
D_G = @(eps_p,tau,s,P,T) eps_p/tau^2*(1-s).^3.*(T/T_ref).^1.5*(P_ref/P); % [-] scaling factor for gas diffusivities
D_L = @(eps_p,tau,s) eps_p/tau^2*s.^3; % [-] scaling factor for liquid species diffusivities
D_H2O_CO2   = @(eps_p,tau,s,T) D_H2O_CO2_ref*D_G(eps_p,tau,s,P,T); % [m^2/s] effective H2O-CO2 gas phase binary diffusion coefficient
D_H2O_CO    = @(eps_p,tau,s,T) D_H2O_CO_ref*D_G(eps_p,tau,s,P,T); % [m^2/s] effective H2O-CO gas phase binary diffusion coefficient
D_H2O_H2    = @(eps_p,tau,s,T) D_H2O_H2_ref*D_G(eps_p,tau,s,P,T); % [m^2/s] effective H2O-H2 gas phase binary diffusion coefficient
D_CO_CO2    = @(eps_p,tau,s,T) D_CO_CO2_ref*D_G(eps_p,tau,s,P,T); % [m^2/s] effective CO-CO2 gas phase binary diffusion coefficient
D_H2_CO     = @(eps_p,tau,s,T) D_H2_CO_ref*D_G(eps_p,tau,s,P,T); % [m^2/s] effective CO-H2 gas phase binary diffusion coefficient
D_H2_CO2    = @(eps_p,tau,s,T) D_H2_CO2_ref*D_G(eps_p,tau,s,P,T); % [m^2/s] effective CO2-H2 gas phase binary diffusion coefficient
D_CO2_L_eff = @(eps_p,tau,s,T) D_CO2_L(T).*D_L(eps_p,tau,s); % [m^2/s] effective CO2 liquid phase binary diffusion coefficient
D_OH_eff    = @(eps_p,tau,s,T) D_OH(T).*D_L(eps_p,tau,s); % [m^2/s] effective OH^- liquid phase binary diffusion coefficient
D_HCO3_eff  = @(eps_p,tau,s,T) D_HCO3(T).*D_L(eps_p,tau,s); % [m^2/s] effective HCO3^- liquid phase binary diffusion coefficient
D_CO3_eff   = @(eps_p,tau,s,T) D_CO3(T).*D_L(eps_p,tau,s); % [m^2/s] effective CO3^2- liquid phase binary diffusion coefficient
x_H2O_GC = RH*P_sat_o(T_C)/P; % [-] mole fraction of water vapor in gas channel
x_CO2_GC = alpha_CO2*(1-x_H2O_GC); % [-] mole fraction of carbon dioxide in gas channel
x_CO_GC = alpha_CO*(1-x_H2O_GC); % [-] mole fraction of carbon monoxide in gas channel
x_H2_GC = alpha_H2*(1-x_H2O_GC); % [-] mole fraction of hydrogen in gas channel
C_CO2_L_b = H_CO2(T_L)*P*x_CO2_GC; % [mol/m^3]
[C_OH_b,C_HCO3_b,C_CO3_b] = carbonation(C_I,C_CO2_L_b,Kw(T_C),K3(T_C),K4(T_C),T_C); % [mol/m^3] dissolved anion species concentrations in bulk electrolyte


% AUXILIARY FUNCTIONS
iff = @(cond,a,b) cond.*a + ~cond.*b; % vectorized ternary operator

% MATERIAL CONSTITUTIVE RELATIONSHIPS
load('GDE_PC_(GDL-Toray)(MPL)(CL)','GDE')
% S_PC = @(P_C,layer,theta) interp2(GDE.(layer).PC , GDE.(layer).theta, GDE.(layer).S , P_C, theta);
% kappa_L_eff = @(kappa,P_C,layer,theta) kappa*(1e-5+interp2(GDE.(layer).PC, GDE.(layer).theta, GDE.(layer).kappa_r_L, P_C, theta)); 
S_PC = @(P_C,layer,theta) interp1(GDE.(layer).PC, GDE.(layer).S(1,:) , P_C);
kappa_L_eff = @(kappa,P_C,layer,theta) kappa*(1e-5+interp1(GDE.(layer).PC, GDE.(layer).kappa_r_L(1,:), P_C)); 
r_K = @(P_C,layer,theta) (1e-6+interp2(GDE.(layer).PC, GDE.(layer).theta, GDE.(layer).r_K, P_C, theta)); % [m] Radius for Knudsen diffusion
s_red = @(s,s_im) (s-s_im)/(1-s_im); % reduced liquid water saturation
gamma_ec = @(x_H2O,x_sat,s,s_im,T) 2e6*iff(x_H2O<x_sat,5e-4*s_red(s,s_im),6e-3*(1-s_red(s,s_im))).*sqrt(R*T/(2*pi*M_w)); % [1/s] evaporation/condensation rate

% INITIAL MESH
Lsum = [0 cumsum(L)];
Nd = numel(L); % number of domains
x = interp1(0:Nd, Lsum, linspace(0, Nd, Nd*10+1));
x = sort([x Lsum(2:end-1)]); % duplicate interface nodes

% SOLVER PREPARATION
sol = bvpinit(x, @yinit);
options = bvpset('Vectorized', 'on', 'NMax', 1e3, 'RelTol', 1e-4, 'AbsTol', 1e-6);

% PARAMETER SWEEP
SOL = cell(size(U));
Np = numel(U); % number of parameters in the sweep
Neq = size(sol.y,1)/2; % number of 2nd-order differential equations
for k = 1:Np
    sol = bvp4c(@odefun, @(ya,yb) bcfun(ya, yb, U(k)), sol, options);
    SOL{k} = sol;
    I(k) = sol.y(14,1)/10; % current density in [mA/cm^2]
end
UI = [U(:) I(:)];
% POSTPROCESSING
Nref = 2; % number of refinements for smoother curve plotting
domains = [1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1 ];
shift = 1e-10;
for k = 1:Np
    x = [];
    for m = 1:Nd
        xa = find(SOL{k}.x==Lsum(m  ), 1, 'last' );
        xb = find(SOL{k}.x==Lsum(m+1), 1, 'first');
        N = xb-xa;
        
        % grid refinement
        x = [x interp1(linspace(0,1,N+1), SOL{k}.x(xa:xb), linspace(shift, 1-shift, N*2^Nref+1))];
        
        % fill solution on inactive domains with NaN
        SOL{k}.y(~kron(domains(:,m),ones(2,1)),xa:xb) = NaN;
    end
    [SOL{k}.y, SOL{k}.yp] = deval(SOL{k}, x);
    SOL{k}.x = x;
end

% PLOT SOLUTION
fig_names = {'Potentials', 'Fluxes'};
unit_scale = [1e-3 1 1 1 1 1 1   1 1 1 1;
              1    1 1 1 1 1 0.1 1 1 1 1 ];
quantity = {'{\itp}_L [kPa]','{\itx}_{H_2O}','{\itx}_{CO_2}','{\itx}_{CO}','{\itx}_{H_2}','{\itT}','{\it\phi}_e [V]','{\itC}_{CO2}','{\itC}_{OH}','{\itC}_{HCO3}','{\itC}_{CO3}';
            '{\itj}_L','{\itj}_{H2O}','{\itj}_{CO_2}','{\itj}_{CO}','{\itj}_{H_2}','{\itj}_{T}','{\itj}_e [mA/cm^2]','{\itj}_{CO_2,L}','{\itj}_{OH}','{\itj}_{HCO3}','{\itj}_{CO3}'};
c = winter(Np);
for m = 1:2
    figure('Name', fig_names{m})
    for n = 1:Neq
        subplot(4,3,n)
        box on
        hold on
        us = unit_scale(m,n);
        for k = 1:Np
            plot(SOL{k}.x*1e6, SOL{k}.y(2*(n-1)+m,:)*us, 'Color', c(k,:), 'DisplayName', [num2str(U(k)),' V'])
        end
        xlim([Lsum(find(domains(n,:),1,'first')) Lsum(find(domains(n,:),1,'last')+1)]*1e6)
        ylim(ylim)
        xlabel('x [μm]')
        ylabel(quantity(m,n))
        for x = Lsum(2:end-1)
            l = line([x x]*1e6, ylim, 'Color', 'k');
            set(get(get(l, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
        end
    end
end

% PLOT POLARIZATION CURVE
figure('Name', 'Polarization curve')
fnplt(cscvn([U; I]))
xlabel({'Cell voltage [V]'})
ylabel('Current density [mA/cm^2]')
set(gca,'XDir','reverse');

% PLOT CAPILLARY PRESSURE PROFILE
figure('Name','Capillary Pressure')
box on
hold on
for k = 1:Np
    plot(SOL{k}.x(1,:)*1e6,SOL{k}.y(1,:)-P,'Color', c(k,:));
end
xlabel('{\itx} [μm]')
ylabel('{\itP}_C (Pa)')

% PLOT SATURATION PROFILE
figure('Name','Saturation')
box on
hold on
for k = 1:Np
    [~,ind] = min(abs(SOL{k}.x(1,:)-Lsum(2)));
    plot([SOL{k}.x(1,1:ind),SOL{k}.x(1,ind+1:end)]*1e6,...
        [S_PC(SOL{k}.y(1,1:ind)-P,'GDL',theta_GDL),S_PC(SOL{k}.y(1,ind+1:end)-P,'MPL',theta_MPL)],'Color', c(k,:));
end
xlabel('{\itx} [μm]')
ylabel('Saturation (-)')

function dydx = odefun(x, y, subdomain)

% READ POTENTIALS & FLUXES
p_L     = y( 1,:); j_L     = y( 2,:);
x_w     = y( 3,:); j_x_w   = y( 4,:);
x_CO2   = y( 5,:); j_x_CO2 = y( 6,:);
x_CO    = y( 7,:); j_x_CO  = y( 8,:);
x_H2    = y( 9,:); j_x_H2  = y(10,:);
T       = y(11,:); j_T     = y(12,:);
phi_e   = y(13,:); j_e     = y(14,:);
C_CO2_L = y(15,:); j_CO2_L = y(16,:);
C_OH    = y(17,:); j_OH    = y(18,:);
C_HCO3  = y(19,:); j_HCO3  = y(20,:);
C_CO3   = y(21,:); j_CO3   = y(22,:);

% ZERO-INITIALIZE ALL DERIVATIVES
z = zeros(size(x));
dp_L     = z; dj_L     = z;
dx_w     = z; dj_x_w   = z;
dx_CO2   = z; dj_x_CO2 = z;
dx_CO    = z; dj_x_CO  = z;
dx_H2    = z; dj_x_H2  = z;
dT       = z; dj_T     = z;
dphi_e   = z; dj_e     = z;
dC_CO2_L = z; dj_CO2_L = z;
dC_OH    = z; dj_OH    = z;
dC_HCO3  = z; dj_HCO3  = z;
dC_CO3   = z; dj_CO3 = z;

% COMPUTE DERIVATIVES
switch subdomain
    case 1 % GAS DIFFUSION LAYER
        p_G = P; % [Pa] gas phase absolute pressure (isobaric)
        p_C = p_L - p_G; % [Pa] capillary pressure
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T,p_C)./P; % saturation water vapor mole fraction
        s = S_PC(p_C,'GDL',theta_GDL); % [-] liquid saturation fraction of pore space
        S_ec = gamma_ec(x_w,x_sat,s,s_im_GDL,T).*C.*(x_w-x_sat); % water evaporation/condensation reaction rate
        S_T = -j_e.*dphi_e + H_ec*S_ec; % Joule heating & evaporation-condensation enthalpy
        dp_L = -j_L./((rho_w/M_w).*kappa_L_eff(kappa_L_GDL,p_C,'GDL',theta_GDL)./mu(T)); % liquid water flux: j_L = -rho_L*(kappa_L/mu_L)*grad(p_L)
        dx_w = d_x_SM(C,x_w,j_x_w,[x_CO2;x_CO;x_H2],[j_x_CO2;j_x_CO;j_x_H2],[D_H2O_CO2(eps_p_GDL,tau_GDL,s,T);D_H2O_CO(eps_p_GDL,tau_GDL,s,T);D_H2O_H2(eps_p_GDL,tau_GDL,s,T)]);     
        dx_CO2 = d_x_SM(C,x_CO2,j_x_CO2,[x_w;x_CO;x_H2],[j_x_w;j_x_CO;j_x_H2],[D_H2O_CO2(eps_p_GDL,tau_GDL,s,T);D_CO_CO2(eps_p_GDL,tau_GDL,s,T);D_H2_CO2(eps_p_GDL,tau_GDL,s,T)]);     
        dx_CO = d_x_SM(C,x_CO,j_x_CO,[x_w;x_CO2;x_H2],[j_x_w;j_x_CO2;j_x_H2],[D_H2O_CO(eps_p_GDL,tau_GDL,s,T);D_CO_CO2(eps_p_GDL,tau_GDL,s,T);D_H2_CO(eps_p_GDL,tau_GDL,s,T)]);
        dx_H2 = d_x_SM(C,x_H2,j_x_H2,[x_w;x_CO2;x_CO],[j_x_w;j_x_CO2;j_x_CO],[D_H2O_H2(eps_p_GDL,tau_GDL,s,T);D_H2_CO2(eps_p_GDL,tau_GDL,s,T);D_H2_CO(eps_p_GDL,tau_GDL,s,T)]);
        dT = -j_T/k_GDL; % heat flux: j_T = -k*grad(T)
        dphi_e = -j_e/sigma_e_GDL; % electron flux: j_e = -sigma_e*grad(phi_e)
        dj_L = S_ec; % conservation of liquid water: div(j_L) = S_L
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_g) = S_H2O
        dj_T = S_T; % conservation of heat: div(j_T) = S_T
    case 2 % MICROPOROUS LAYER
        p_G = P; % [Pa] gas phase absolute pressure (isobaric)
        p_C = p_L - p_G; % [Pa] capillary pressure
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T,p_C)./P; % saturation water vapor mole fraction
        s = S_PC(p_C,'MPL',theta_MPL); % [-] liquid saturation fraction of pore space
        S_ec = gamma_ec(x_w,x_sat,s,s_im_MPL,T).*C.*(x_w-x_sat); % water evaporation/condensation reaction rate
        S_T = -j_e.*dphi_e + H_ec*S_ec; % Joule heating & evaporation-condensation enthalpy
        dp_L = -j_L./((rho_w/M_w).*kappa_L_eff(kappa_L_MPL,p_C,'MPL',theta_MPL)./mu(T)); % liquid water flux: j_L = -rho_L*(kappa_L/mu_L)*grad(p_L)
        dx_w = d_x_SM(C,x_w,j_x_w,[x_CO2;x_CO;x_H2],[j_x_CO2;j_x_CO;j_x_H2],[D_H2O_CO2(eps_p_MPL,tau_MPL,s,T);D_H2O_CO(eps_p_MPL,tau_MPL,s,T);D_H2O_H2(eps_p_MPL,tau_MPL,s,T)]);     
        dx_CO2 = d_x_SM(C,x_CO2,j_x_CO2,[x_w;x_CO;x_H2],[j_x_w;j_x_CO;j_x_H2],[D_H2O_CO2(eps_p_MPL,tau_MPL,s,T);D_CO_CO2(eps_p_MPL,tau_MPL,s,T);D_H2_CO2(eps_p_MPL,tau_MPL,s,T)]);     
        dx_CO = d_x_SM(C,x_CO,j_x_CO,[x_w;x_CO2;x_H2],[j_x_w;j_x_CO2;j_x_H2],[D_H2O_CO(eps_p_MPL,tau_MPL,s,T);D_CO_CO2(eps_p_MPL,tau_MPL,s,T);D_H2_CO(eps_p_MPL,tau_MPL,s,T)]);
        dx_H2 = d_x_SM(C,x_H2,j_x_H2,[x_w;x_CO2;x_CO],[j_x_w;j_x_CO2;j_x_CO],[D_H2O_H2(eps_p_MPL,tau_MPL,s,T);D_H2_CO2(eps_p_MPL,tau_MPL,s,T);D_H2_CO(eps_p_MPL,tau_MPL,s,T)]);
        dT = -j_T/k_MPL; % heat flux: j_T = -k*grad(T)
        dphi_e = -j_e/sigma_e_MPL; % electron flux: j_e = -sigma_e*grad(phi_e)
        dj_L = S_ec; % conservation of liquid water: div(j_L) = S_L
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_g) = S_H2O
        dj_T = S_T; % conservation of heat: div(j_T) = S_T
    case 3 % CATALYST LAYER
        p_G = P; % [Pa] gas phase absolute pressure (isobaric)
        p_C = p_L - p_G; % [Pa] capillary pressure
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T,p_C)./P; % saturation water vapor mole fraction
        s = S_PC(p_C,'CL',theta_CL); % [-] liquid saturation fraction of pore space
        S_ec = gamma_ec(x_w,x_sat,s,s_im_CL,T).*C.*(x_w-x_sat); % evaporation/condensation reaction rate
        
        k_GL = 0.1; % [m/s] 
        S_pt_CO2 = 2e6*k_GL*(H_CO2(T).*p_G.*x_CO2 - C_CO2_L); % homogeneous mass transfer of carbon dioxide between gas and liquid phases
        S_h_3_f = k3f*C_CO2_L.*C_OH; % ????? do these need to be saturation corrected?
        S_h_4_f = k4f*C_HCO3.*C_OH; % ????? do these need to be saturation corrected?
        S_h_3_r = k3r(T).*C_HCO3; % REVERSE REACTION = COULD CAUSE ISSUES ????? do these need to be saturation corrected?
        S_h_4_r = k4r(T).*C_CO3; % REVERSE REACTION = COULD CAUSE ISSUES ????? do these need to be saturation corrected?
        pH = 14-(-log10(C_OH)); % electrolyte pH
        
        eta_COER = iff((phi_e-U_0_COER)>0,0,phi_e-U_0_COER); % COER overpotential
        eta_HER = iff((phi_e-U_0_HER)>0,0,phi_e-U_0_HER); % HER overpotential
        i_COER = BV(i_0_COER(T),C_CO2_L,C_ref_CO2,gamma_BV_CO2,beta_COER,T,eta_COER); % COER electrochemical reaction rate
        i_HER = BV(i_0_HER(T,pH),p_ref_H2,p_ref_H2,gamma_BV_H,beta_HER,T,eta_HER); % HER electrochemical reaction rate
        S_H2O = -a_CCL*(i_COER/(2*F)+i_HER/F); % liquid water reaction rate (Faraday's law)
        S_CO2 = -S_pt_CO2; % volumetric dissolution of carbon dioxide in liquid phase
        S_CO = a_CCL*i_COER/(2*F); % carbon monoxide reaction rate (Faraday's law)
        S_H2 = a_CCL*i_HER/(2*F); % hydrogen reaction rate (Faraday's law)
        S_T = -j_e.*dphi_e -a_CCL*i_COER.*(eta_COER+Pi_COER(T))-a_CCL*i_HER.*(eta_HER+Pi_HER(T)) + H_ec*S_ec; % Joule+Peltier heating & evaporation-condensation enthalpy
        S_e = -a_CCL*(i_COER+i_HER);
        
        S_CO2_L = S_pt_CO2 - S_h_3_f - a_CCL*i_COER/(2*F); % carbon dioxide absorption from gas phase and electrochemical reaction rate (Faraday's law)
        S_OH = a_CCL*(i_COER/F+i_HER/F) - S_h_3_f - 1e-7*S_h_4_f;
        S_HCO3 = (S_h_3_f - 1e-7*S_h_4_f) + (S_h_4_r - S_h_3_r);
        S_CO3 = S_h_4_f - S_h_4_r;
        
%         S_h_3_f - S_h_3_r
%         S_h_4_f - S_h_4_r
        
        dp_L = -j_L./((rho_w/M_w).*kappa_L_eff(kappa_L_CL,p_C,'CL',theta_CL)./mu(T)); % liquid water flux: j_L = -rho_L*(kappa_L/mu_L)*grad(p_L)
        dx_w = d_x_SM(C,x_w,j_x_w,[x_CO2;x_CO;x_H2],[j_x_CO2;j_x_CO;j_x_H2],[D_H2O_CO2(eps_p_CL,tau_CL,s,T);D_H2O_CO(eps_p_CL,tau_CL,s,T);D_H2O_H2(eps_p_CL,tau_CL,s,T)]);     
        dx_CO2 = d_x_SM(C,x_CO2,j_x_CO2,[x_w;x_CO;x_H2],[j_x_w;j_x_CO;j_x_H2],[D_H2O_CO2(eps_p_CL,tau_CL,s,T);D_CO_CO2(eps_p_CL,tau_CL,s,T);D_H2_CO2(eps_p_CL,tau_CL,s,T)]);     
        dx_CO = d_x_SM(C,x_CO,j_x_CO,[x_w;x_CO2;x_H2],[j_x_w;j_x_CO2;j_x_H2],[D_H2O_CO(eps_p_CL,tau_CL,s,T);D_CO_CO2(eps_p_CL,tau_CL,s,T);D_H2_CO(eps_p_CL,tau_CL,s,T)]);
        dx_H2 = d_x_SM(C,x_H2,j_x_H2,[x_w;x_CO2;x_CO],[j_x_w;j_x_CO2;j_x_CO],[D_H2O_H2(eps_p_CL,tau_CL,s,T);D_H2_CO2(eps_p_CL,tau_CL,s,T);D_H2_CO(eps_p_CL,tau_CL,s,T)]);
        dT = -j_T/k_CL; % heat flux: j_T = -k*grad(T)
        dphi_e = -j_e/sigma_e_CL; % electron flux: j_e = -sigma_e*grad(phi_e)
        
        dC_CO2_L = -j_CO2_L./D_CO2_L_eff(eps_p_CL,tau_CL,s,T); % dissolved CO2 flux: j_CO2_L = -D*grad(C_CO2_L)
        dC_OH = -j_OH./D_OH_eff(eps_p_CL,tau_CL,s,T); % dissolved OH flux: j_OH = -D*grad(C_OH)
        dC_HCO3 = -j_HCO3./D_HCO3_eff(eps_p_CL,tau_CL,s,T); % dissolved HCO3 flux: j_HCO3 = -D*grad(C_HCO3)
        dC_CO3 = -j_CO3./D_CO3_eff(eps_p_CL,tau_CL,s,T); % dissolved CO3 flux: j_CO3 = -D*grad(C_CO3)
        
        dj_L = S_ec + S_H2O; % conservation of liquid water: div(j_L) = S_L
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_v) = S_H2O
        dj_x_CO2 = S_CO2; % conservation of carbon dioxide gas: div(j_CO2) = S_CO2
        dj_x_CO = S_CO; % conservation of carbon monoxide gas: div(j_CO) = S_CO
        dj_x_H2 = S_H2; % conservation of hydrogen gas: div(j_CO) = S_CO
        dj_T = S_T; % conservation of heat: div(j_T) = S_T
        dj_e = S_e; % conservation of electrons: div(j_e) = S_e
        
        dj_CO2_L = S_CO2_L; % conservation of dissolved carbon dioxide: div(j_CO2_L) = S_CO2_L
        dj_OH = S_OH; % conservation of hydroxide anions: div(j_OH) = S_OH
        dj_HCO3 = S_HCO3; % conservation of bicarbonate anions: div(j_CO2_L) = S_HCO3
        dj_CO3 = S_CO3; % conservation of carbonate dianions: div(j_CO3) = S_CO3
end

% ASSEMBLE DERIVATIVES
dydx = [dp_L     ; dj_L     ;
        dx_w     ; dj_x_w   ;
        dx_CO2   ; dj_x_CO2 ;
        dx_CO    ; dj_x_CO  ;
        dx_H2    ; dj_x_H2  ;
        dT       ; dj_T     ;
        dphi_e   ; dj_e     ;
        dC_CO2_L ; dj_CO2_L ;
        dC_OH    ; dj_OH    ;
        dC_HCO3  ; dj_HCO3  ;
        dC_CO3   ; dj_CO3    ];
end

function y0 = yinit(x, subdomain)
% POTENTIALS INITIALLY SET TO GAS CHANNEL CONDITIONS
p_L = p_C_GC + P;
x_w = x_H2O_GC;
x_CO2 = x_CO2_GC;
x_CO = x_CO_GC;
x_H2 = x_H2_GC;
T = T_C;
phi_e = U(1);
C_CO2_L = C_CO2_L_b;
C_OH = C_OH_b;
C_HCO3 = C_HCO3_b;
C_CO3 = C_CO3_b;

% ALL FLUXES ARE INITIALLY ZERO
y0 = [p_L     ; 0 ; 
      x_w     ; 0 ;
      x_CO2   ; 0 ;
      x_CO    ; 0 ;
      x_H2    ; 0 ;      
      T       ; 0 ;
      phi_e   ; 0 ;
      C_CO2_L ; 0 ;
      C_OH    ; 0 ;
      C_HCO3  ; 0 ;
      C_CO3   ; 0  ];
end

function res = bcfun(ya, yb, U)

res = ya(:); % homogeneous BC everywhere by default

% LIQUID WATER

res(0*Neq+1) = ya(1,1) - (p_C_GC + P); % GC liquid pressure
res(0*Neq+2) = yb(2,1) - ya(2,2); % flux continuity between GDL & MPL  
res(2*Neq+1) = ya(1,2) - yb(1,1); % potential continuity between GDL & MPL
res(2*Neq+2) = yb(2,2) - ya(2,3); % flux continuity between MPL & CL 
res(4*Neq+1) = yb(1,2) - ya(1,3); % potential continuity between MPL & CL
res(4*Neq+2) = yb(1,3) - (p_C_LC + P); % LC liquid pressure

% WATER VAPOR
res(0*Neq+3) = ya(3,1) - x_H2O_GC; % gas channel water vapor content
res(0*Neq+4) = yb(4,1) - ya(4,2); % flux continuity between GDL & MPL
res(2*Neq+3) = yb(3,1) - ya(3,2); % potential continuity between GDL & MPL
res(2*Neq+4) = yb(4,2) - ya(4,3); % flux continuity between MPL & CL
res(4*Neq+3) = yb(3,2) - ya(3,3); % potential continuity between MPL & CL
res(4*Neq+4) = yb(4,3); % zero flux at boundary
% res(4*Neq+4) = yb(3,3) - P_sat_o(T_L)/P; % saturated water vapor at liquid channel

% CARBON DIOXIDE GAS
res(0*Neq+5) = ya(5,1) - x_CO2_GC; % gas channel carbon dioxide content
res(0*Neq+6) = yb(6,1) - ya(6,2); % flux continuity between GDL & MPL
res(2*Neq+5) = yb(5,1) - ya(5,2); % potential continuity between GDL & MPL
res(2*Neq+6) = yb(6,2) - ya(6,3); % flux continuity between MPL & CL
res(4*Neq+5) = yb(5,2) - ya(5,3); % potential continuity between MPL & CL
res(4*Neq+6) = yb(6,3); % zero flux at boundary

% CARBON MONOXIDE GAS
res(0*Neq+7) = ya(7,1) - x_CO_GC; % gas channel carbon monoxide content
res(0*Neq+8) = yb(8,1) - ya(8,2); % flux continuity between GDL & MPL
res(2*Neq+7) = yb(7,1) - ya(7,2); % potential continuity between GDL & MPL
res(2*Neq+8) = yb(8,2) - ya(8,3); % flux continuity between MPL & CL
res(4*Neq+7) = yb(7,2) - ya(7,3); % potential continuity between MPL & CL
res(4*Neq+8) = yb(8,3); % zero flux at boundary

% HYDROGEN GAS
res(0*Neq+ 9) = ya( 9,1) - x_H2_GC; % gas channel hydrogen content
res(0*Neq+10) = yb(10,1) - ya(10,2); % flux continuity between GDL & MPL
res(2*Neq+ 9) = yb( 9,1) - ya( 9,2); % potential continuity between GDL & MPL
res(2*Neq+10) = yb(10,2) - ya(10,3); % flux continuity between MPL & CL
res(4*Neq+ 9) = yb( 9,2) - ya( 9,3); % potential continuity between MPL & CL
res(4*Neq+10) = yb(10,3); % zero flux at boundary

% TEMPERATURE
res(0*Neq+11) = ya(11,1) - T_C; % cathode boundary temperature
res(0*Neq+12) = yb(11,3) - T_L; % liquid boundary temperature
for d = 2:Nd
    res(2*(d-1)*Neq+11) = ya(11,d) - yb(11,d-1); % potential continuity
    res(2*(d-1)*Neq+12) = ya(12,d) - yb(12,d-1); % flux continuity
end

% ELECTRONS
res(0*Neq+13) = ya(13,1) - U; % current collector electrical potential
res(0*Neq+14) = yb(14,1) - ya(14,2); % flux continuity between GDL & MPL
res(2*Neq+13) = yb(13,1) - ya(13,2); % potential continuity between GDL & MPL
res(2*Neq+14) = yb(14,2) - ya(14,3); % flux continuity between MPL & CL
res(4*Neq+13) = yb(13,2) - ya(13,3); % potential continuity between MPL & CL
res(4*Neq+14) = yb(14,3); % zero flux at boundary
    
% CARBON DIOXIDE DISSOLVED IN LIQUID
res(0*Neq+15) = ya(15,1) - H_CO2(T_C)*P*x_CO2_GC; % Dissolved CO2 concentration at GC interface
res(0*Neq+16) = yb(16,1) - ya(16,2); % flux continuity between GDL & MPL
res(2*Neq+15) = yb(15,1) - ya(15,2); % potential continuity between GDL & MPL
res(2*Neq+16) = yb(16,2) - ya(16,3); % flux continuity between MPL & CL
res(4*Neq+15) = yb(15,2) - ya(15,3); % potential continuity between MPL & CL
res(4*Neq+16) = yb(15,3) - C_CO2_L_b; % Dissolved CO2 concentration in bulk electrolyte

% HYDROXIDE ANIONS
res(0*Neq+17) = ya(18,1); % no flux at boundary
res(0*Neq+18) = yb(18,1) - ya(18,2); % flux continuity between GDL & MPL
res(2*Neq+17) = yb(17,1) - ya(17,2); % potential continuity between GDL & MPL
res(2*Neq+18) = yb(18,2) - ya(18,3); % flux continuity between MPL & CL
res(4*Neq+17) = yb(17,2) - ya(17,3); % potential continuity between MPL & CL
res(4*Neq+18) = yb(17,3) - C_OH_b; % Concentration in bulk electrolyte

% BICARBONATE ANIONS
res(0*Neq+19) = ya(20,1); % no flux at boundary
res(0*Neq+20) = yb(20,1) - ya(20,2); % flux continuity between GDL & MPL
res(2*Neq+19) = yb(19,1) - ya(19,2); % potential continuity between GDL & MPL
res(2*Neq+20) = yb(20,2) - ya(20,3); % flux continuity between MPL & CL
res(4*Neq+19) = yb(19,2) - ya(19,3); % potential continuity between MPL & CL
res(4*Neq+20) = yb(19,3) - C_HCO3_b; % Concentration in bulk electrolyte

% CARBONATE DI-ANIONS
res(0*Neq+21) = ya(22,1); % no flux at boundary
res(0*Neq+22) = yb(22,1) - ya(22,2); % flux continuity between GDL & MPL
res(2*Neq+21) = yb(21,1) - ya(21,2); % potential continuity between GDL & MPL
res(2*Neq+22) = yb(22,2) - ya(22,3); % flux continuity between MPL & CL
res(4*Neq+21) = yb(21,2) - ya(21,3); % potential continuity between MPL & CL
res(4*Neq+22) = yb(21,3) - C_CO3_b; % Concentration in bulk electrolyte

end

function d_x_i = d_x_SM(C_T,x_i,N_i,X_j,N_j,D_ij)
    % m   : number of domain position elements [-] (1-D domain vector length)
    % C_T : total gas concentration scalar [mol/m^3]
    % x_i : mole fraction vector component i [-], (1 x m)
    % N_i : molar flux vector component i [-], (1 x m)
    % X_j : mole fraction matrix components j [-], (1+count(j)) x m)
    % N_j : molar flux matrix components j [-], (1+count(j)) x m)
    % D_ij: binary diffusivity vectors for {i,j} [m^2/s], (count(j) x m)
    d_x_i = zeros(1,length(x_i));
    for j = 1:size(X_j,1)
        d_x_i = d_x_i + (x_i .* N_j(j,:) - X_j(j,:) .* N_i) ./ (C_T .* D_ij(j,:));
    end
end

function [C_OH,C_HCO3,C_CO3] = carbonation(C_cation,C_CO2_AQ,Kw,K1,K2,T)
    % C_cation: [mol/m^3] cation concentration
    % p_CO2: [Pa] partial pressure of carbon dioxide gas
    % T: [K] temperature
    % C_OH: [mol/m^3] equilibrium hydroxide concentration
    % C_HCO3: [mol/m^3] equilibrium bicarbonate concentration
    % C_CO3: [mol/m^3] equilibrium carbonate concentration
    fun = @(C) objective(C,K1,K2,Kw,C_CO2_AQ,C_cation);
    C_0 = [ 0.0001 ; C_cation ; 0.0001 ];
    C_eqm = fsolve(fun,C_0);
    C_OH = C_eqm(1,:);
    C_HCO3 = C_eqm(2,:);
    C_CO3 = C_eqm(3,:);
    function F = objective(C,K1,K2,Kw,C_CO2,C_cation)
        C_OH = C(1);
        C_HCO3 = C(2);
        C_CO3 = C(3);
        F(1) = K1 * C_CO2 * C_OH - C_HCO3 ;
        F(2) = K2 * C_HCO3 * C_OH - C_CO3 ;
        F(3) = C_cation + ( Kw / C_OH - C_OH ) - C_HCO3 - 2 * C_CO3 ; 
    end
end

end