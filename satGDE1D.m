function [sol] = satGDE1D

% CONSTANTS
M_w = 18e-3; % [kg/mol] molar mass of water
rho_w = 997; % [kg/m^3] density of water at reference temperature
T_0 = 273.15; % [K] zero degrees Celsius
R = 8.31446; % [J*mol/K] universal gas constant
F = 96485.333; % [C/mol] Faraday constant
P_ref = 101325; % [Pa] reference pressure
T_ref = T_0+25; % [K] reference temperature

% OPERATING CONDITIONS
P = 1.5e5; % [Pa] total pressure in gas channel
RH = 1.0; % [-] relative humidity in gas channel
p_C_GC = 0; % [Pa] capillary pressure at GDL/GC interface
p_C_LC = 500; % [Pa] capillary pressure at Liquid channel interface
T_C = T_0+30; % [K] temperature of gas channel
T_L = T_0+30; % [K] temperature of liquid channel
alpha_CO2 = 0.9999; % [-] mole fraction of carbon dioxide in dry feed gas
alpha_CO = 0.0001; % [-] mole fraction of carbon monoxide in dry feed gas

N_w = 0.0; % [mol/(m^2*s)] water injection flux through CL
i = [0:2:20]'; % [A/m^2] current density

% MATERIAL PARAMETERS
% L = [300 45]*1e-6; % [m] gas diffusion electrode domain thicccnesses
L = [300 45 10]*1e-6; % [m] gas diffusion electrode domain thicccnesses
a_CCL = 1e7; % [m^2/m^3] ECSA density of CCL
H_ec = 42e3; % [J/mol] molar enthalphy of evaporation/condensation
k_GDL = 1.2; % [W/(m*K)] thermal conductivity of GDL
k_MPL = 0.2; % [W/(m*K)] thermal conductivity of MPL
k_CL = 0.27; % [W/(m*K)] thermal conductivity of CL
s_im_GDL = 0.05; % [-] immobile liquid water saturation of GDL
s_im_MPL = 1e-9; % [-] immobile liquid water saturation of MPL
s_im_CL = 0.05; % [-] immobile liquid water saturation of CL
eps_p_GDL = 0.7; % [-] porosity of GDL
eps_p_MPL = 0.3; % [-] porosity of MPL
eps_p_CL = 0.4; % [-] porosity of CL
kappa_L_GDL = 0.8e-11; % [m^2] absolute permeability of GDL
kappa_L_MPL = 5e-14; % [m^2] absolute permeability of MPL
kappa_L_CL = 1e-13; % [m^2] absolute permeability of CL
tau_GDL = 1.6; % [-] pore tortuosity of GDL
tau_MPL = 1.6; % [-] pore tortuosity of MPL
tau_CL = 1.6; % [-] pore tortuosity of CL
theta_GDL = 93; % [°] intrinsic mean contact angle of GDL
theta_MPL = 110; % [°] intrinsic mean contact angle of MPL
theta_CL = 93; % [°] intrinsic mean contact angle of CL

% WATER CONSTITUTIVE RELATIONSHIPS
P_sat_o = @(T) exp(23.1963-3816.44./(T-46.13)); % [Pa] uncorrected saturation pressure of water vapor
P_sat = @(T,P_C) P_sat_o(T).*exp(P_C.*(M_w/rho_w)./(R.*T)); % [Pa] capillary pressure corrected saturation pressure of water vapor
mu = @(T) 1e-3*exp(-3.63148+542.05./(T-144.15)); % [Pa*s] dynamic viscosity of liquid water

% DIFFUSION COEFFICIENTS (BULK)
D_ab = @(nu_p_a,nu_p_b,M_a,M_b,P,T) ...
    1e-4*(1e-3*T.^1.75.*(1/M_a+1/M_b)^0.5)./(P/101325*(nu_p_a^0.33+nu_p_b^0.33)^2); %[m^2/s] binary diffusion coefficient
nu_p_H2O = 12.7; % [-] diffusion volume of gaseous species
nu_p_CO2 = 26.9; % [-] diffusion volume of gaseous species
nu_p_CO = 18.9; % [-] diffusion volume of gaseous species
M_H2O = 18.02; % [g/mol] molar mass of gaseous species
M_CO2 = 44.01; % [g/mol] molar mass of gaseous species
M_CO = 28.01; % [g/mol] molar mass of gaseous species
D_H2O_CO2_ref = D_ab(nu_p_H2O,nu_p_CO2,M_H2O,M_CO2,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, H2O in CO2
D_H2O_CO_ref = D_ab(nu_p_H2O,nu_p_CO,M_H2O,M_CO,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, H2O in CO2
D_CO_CO2_ref = D_ab(nu_p_CO,nu_p_CO2,M_CO,M_CO2,P_ref,T_ref); % [m^2/s] reference diffusion coefficient, CO in CO2

% MODEL PARAMETERIZATION
D = @(eps_p,tau,s,P,T) eps_p/tau^2*(1-s).^3.*(T/T_ref).^1.5*(P_ref/P); % [-] scaling factor for gas diffusivities
D_H2O_CO2 = @(eps_p,tau,s,T) D_H2O_CO2_ref*D(eps_p,tau,s,P,T); % [m^2/s] H2O gas phase diffusion coefficient
D_H2O_CO = @(eps_p,tau,s,T) D_H2O_CO_ref*D(eps_p,tau,s,P,T); % [m^2/s] H2O gas phase diffusion coefficient
D_CO_CO2 = @(eps_p,tau,s,T) D_CO_CO2_ref*D(eps_p,tau,s,P,T); % [m^2/s] CO2 gas phase diffusion coefficient
x_H2O_GC = RH*P_sat_o(T_C)/P; % [-] mole fraction of water vapor in gas channel
x_CO2_GC = alpha_CO2*(1-x_H2O_GC); % [-] mole fraction of carbon dioxide in gas channel
x_CO_GC = alpha_CO*(1-x_H2O_GC); % [-] mole fraction of carbon monoxide in gas channel

% AUXILIARY FUNCTIONS
iff = @(cond,a,b) cond.*a + ~cond.*b; % vectorized ternary operator

% MATERIAL CONSTITUTIVE RELATIONSHIPS
% load('S_PC_(GDL-Toray)(MPL)(CL).mat','SatPC');
% load('kappa_PC_(GDL-Toray)(MPL)(CL).mat','kappaPC');
load('GDE_PC_(GDL-Toray)(MPL)(CL)','GDE')
% S_PC = @(P_C,layer) interp1(SatPC.(layer).PC , SatPC.(layer).S , P_C);
% kappa_eff_L = @(kappa,P_C,layer) kappa*(1e-9+interp1(kappaPC.(layer).PC,kappaPC.(layer).kappa_r_L,P_C));
S_PC = @(P_C,layer,theta) interp2(GDE.(layer).PC , GDE.(layer).theta, GDE.(layer).S , P_C, theta);
kappa_L_eff = @(kappa,P_C,layer,theta) kappa*(1e-5+interp2(GDE.(layer).PC, GDE.(layer).theta, GDE.(layer).kappa_r_L, P_C, theta)); 
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
SOL = cell(size(i));
Np = numel(i); % number of parameters in the sweep
Neq = size(sol.y,1)/2; % number of 2nd-order differential equations
for k = 1:Np
    sol = bvp4c(@odefun, @(ya,yb) bcfun(ya, yb, i(k)), sol, options);
    SOL{k} = sol;
end
% POSTPROCESSING
Nref = 2; % number of refinements for smoother curve plotting
domains = [1 1 1;
           1 1 1;
           1 1 1;
           1 1 1;
           1 1 1];
% domains = [1 1;
%            1 1;
%            1 1;
%            1 1;
%            1 1;];
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
unit_scale = [1 1 1 1 1;
              1 1 1 1 1;
              1 1 1 1 1];
% unit_scale = [1 1 1 1 1;
%               1 1 1 1 1];
quantity = {'{\itp}_L','{\itx}_{H_2O}','{\itx}_{CO_2}','{\itx}_{CO}','{\itT}';
            '{\itj}_L','{\itj}_{H2O}','{\itj}_{CO_2}','{\itj}_{CO}','{\itj}_{T}'};
c = winter(Np);
for m = 1:2
    figure('Name', fig_names{m})
    for n = 1:Neq
        subplot(3,2,n)
        box on
        hold on
        us = unit_scale(m,n);
        for k = 1:Np
            plot(SOL{k}.x*1e6, SOL{k}.y(2*(n-1)+m,:)*us, 'Color', c(k,:), 'DisplayName', [num2str(i(k)),' A/m^2'])
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
    legend(strcat(cellstr(num2str(i)),' A/m^2'),'Location','best');
end

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
p_L   = y( 1,:); j_L     = y( 2,:);
x_w   = y( 3,:); j_x_w   = y( 4,:);
x_CO2 = y( 5,:); j_x_CO2 = y( 6,:);
x_CO  = y( 7,:); j_x_CO  = y( 8,:);
T     = y( 9,:); j_T     = y( 10,:);

% ZERO-INITIALIZE ALL DERIVATIVES
z = zeros(size(x));
dp_L   = z; dj_L     = z;
dx_w   = z; dj_x_w   = z;
dx_CO2 = z; dj_x_CO2 = z;
dx_CO  = z; dj_x_CO  = z;
dT     = z; dj_T     = z;

% COMPUTE DERIVATIVES
switch subdomain
    case 1 % GAS DIFFUSION LAYER
        p_G = P; % [Pa] gas phase absolute pressure (isobaric)
        p_C = p_L - p_G; % [Pa] capillary pressure
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T,p_C)./P; % saturation water vapor mole fraction
        s = S_PC(p_C,'GDL',theta_GDL); % [-] liquid saturation fraction of pore space
        S_ec = gamma_ec(x_w,x_sat,s,s_im_GDL,T).*C.*(x_w-x_sat); % water evaporation/condensation reaction rate
        dp_L = -j_L./((rho_w/M_w).*kappa_L_eff(kappa_L_GDL,p_C,'GDL',theta_GDL)./mu(T)); % liquid water flux: j_L = -rho_L*(kappa_L/mu_L)*grad(p_L)
        dx_w = d_x_SM(C,x_w,j_x_w,[x_CO2;x_CO],[j_x_CO2;j_x_CO],[D_H2O_CO2(eps_p_GDL,tau_GDL,s,T);D_H2O_CO(eps_p_GDL,tau_GDL,s,T)]);     
        dx_CO2 = d_x_SM(C,x_CO2,j_x_CO2,[x_w;x_CO],[j_x_w;j_x_CO],[D_H2O_CO2(eps_p_GDL,tau_GDL,s,T);D_CO_CO2(eps_p_GDL,tau_GDL,s,T)]);     
        dx_CO = d_x_SM(C,x_CO,j_x_CO,[x_w;x_CO2],[j_x_w;j_x_CO2],[D_H2O_CO(eps_p_GDL,tau_GDL,s,T);D_CO_CO2(eps_p_GDL,tau_GDL,s,T)]);
        dT = -j_T/k_GDL; % heat flux: j_T = -k*grad(T)
        dj_L = S_ec; % conservation of liquid water: div(j_L) = S_L
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_g) = S_H2O
        dj_T = H_ec*S_ec; % conservation of heat: div(j_T) = S_T
    case 2 % MICROPOROUS LAYER
        p_G = P; % [Pa] gas phase absolute pressure (isobaric)
        p_C = p_L - p_G; % [Pa] capillary pressure
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T,p_C)./P; % saturation water vapor mole fraction
        s = S_PC(p_C,'MPL',theta_MPL); % [-] liquid saturation fraction of pore space
        S_ec = gamma_ec(x_w,x_sat,s,s_im_MPL,T).*C.*(x_w-x_sat); % water evaporation/condensation reaction rate
        dp_L = -j_L./((rho_w/M_w).*kappa_L_eff(kappa_L_MPL,p_C,'MPL',theta_MPL)./mu(T)); % liquid water flux: j_L = -rho_L*(kappa_L/mu_L)*grad(p_L)
        dx_w = d_x_SM(C,x_w,j_x_w,[x_CO2;x_CO],[j_x_CO2;j_x_CO],[D_H2O_CO2(eps_p_MPL,tau_MPL,s,T);D_H2O_CO(eps_p_MPL,tau_MPL,s,T)]);     
        dx_CO2 = d_x_SM(C,x_CO2,j_x_CO2,[x_w;x_CO],[j_x_w;j_x_CO],[D_H2O_CO2(eps_p_MPL,tau_MPL,s,T);D_CO_CO2(eps_p_MPL,tau_MPL,s,T)]);     
        dx_CO = d_x_SM(C,x_CO,j_x_CO,[x_w;x_CO2],[j_x_w;j_x_CO2],[D_H2O_CO(eps_p_MPL,tau_MPL,s,T);D_CO_CO2(eps_p_MPL,tau_MPL,s,T)]);
        dT = -j_T/k_MPL; % heat flux: j_T = -k*grad(T)
        dj_L = S_ec; % conservation of liquid water: div(j_L) = S_L
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_g) = S_H2O
        dj_T = H_ec*S_ec; % conservation of heat: div(j_T) = S_T
    case 3 % CATALYST LAYER
        p_G = P; % [Pa] gas phase absolute pressure (isobaric)
        p_C = p_L - p_G; % [Pa] capillary pressure
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T,p_C)./P; % saturation water vapor mole fraction
        s = S_PC(p_C,'CL',theta_CL); % [-] liquid saturation fraction of pore space
        S_ec = gamma_ec(x_w,x_sat,s,s_im_CL,T).*C.*(x_w-x_sat); % evaporation/condensation reaction rate
        S_H2O = -ones(size(x_w))*a_CCL*i(k)/(2*F); % water reaction rate (Faraday's law)
        S_CO2 = -ones(size(x_CO2))*a_CCL*i(k)/(2*F); % carbon dioxide reaction rate (Faraday's law)
        S_CO = ones(size(x_CO))*a_CCL*i(k)/(2*F); % carbon monoxide reaction rate (Faraday's law)
        dp_L = -j_L./((rho_w/M_w).*kappa_L_eff(kappa_L_CL,p_C,'CL',theta_CL)./mu(T)); % liquid water flux: j_L = -rho_L*(kappa_L/mu_L)*grad(p_L)
        dx_w = d_x_SM(C,x_w,j_x_w,[x_CO2;x_CO],[j_x_CO2;j_x_CO],[D_H2O_CO2(eps_p_CL,tau_CL,s,T);D_H2O_CO(eps_p_CL,tau_CL,s,T)]);     
        dx_CO2 = d_x_SM(C,x_CO2,j_x_CO2,[x_w;x_CO],[j_x_w;j_x_CO],[D_H2O_CO2(eps_p_CL,tau_CL,s,T);D_CO_CO2(eps_p_CL,tau_CL,s,T)]);     
        dx_CO = d_x_SM(C,x_CO,j_x_CO,[x_w;x_CO2],[j_x_w;j_x_CO2],[D_H2O_CO(eps_p_CL,tau_CL,s,T);D_CO_CO2(eps_p_CL,tau_CL,s,T)]);
        dT = -j_T/k_CL; % heat flux: j_T = -k*grad(T)
        dj_L = S_ec + S_H2O; % conservation of liquid water: div(j_L) = S_L
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_v) = S_H2O
        dj_x_CO2 = S_CO2; % conservation of carbon dioxide gas: div(j_CO2) = S_CO2
        dj_x_CO = S_CO; % conservation of carbon monoxide gas: div(j_CO) = S_CO
        dj_T = H_ec*S_ec; % conservation of heat: div(j_T) = S_T
end

% ASSEMBLE DERIVATIVES
dydx = [dp_L   ; dj_L    ;
        dx_w   ; dj_x_w  ;
        dx_CO2 ; dj_x_CO2;
        dx_CO  ; dj_x_CO ;
        dT     ; dj_T    ];
end

function y0 = yinit(x, subdomain)
% POTENTIALS INITIALLY SET TO GAS CHANNEL CONDITIONS
p_L = p_C_GC + P;
x_w = x_H2O_GC;
x_CO2 = x_CO2_GC;
x_CO = x_CO_GC;
T = T_C;

% ALL FLUXES ARE INITIALLY ZERO
y0 = [p_L   ; 0 ; 
      x_w   ; 0 ;
      x_CO2 ; 0 ;
      x_CO  ; 0 ;
      T     ; 0 ;];
end

function res = bcfun(ya, yb, i)

res = ya(:); % homogeneous BC everywhere by default

% LIQUID WATER

res(0*Neq+1) = ya(1,1) - (p_C_GC + P); % GC liquid pressure
res(0*Neq+2) = yb(2,1) - ya(2,2); % flux continuity between GDL & MPL  
res(2*Neq+1) = ya(1,2) - yb(1,1); % potential continuity between GDL & MPL
res(2*Neq+2) = yb(2,2) - ya(2,3); % flux continuity between MPL & CL 
% res(2*Neq+2) = yb(2,2) + i/(2*F); % water injection flux to CL
res(4*Neq+1) = yb(1,2) - ya(1,3); % potential continuity between MPL & CL
% res(4*Neq+2) = yb(2,3) + N_w; % water injection flux to CL
res(4*Neq+2) = yb(1,3) - (p_C_LC + P); % LC liquid pressure

% WATER VAPOR
res(0*Neq+3) = ya(3,1) - x_H2O_GC; % gas channel water vapor content
res(0*Neq+4) = yb(4,1) - ya(4,2); % flux continuity between GDL & MPL
res(2*Neq+3) = yb(3,1) - ya(3,2); % potential continuity between GDL & MPL
res(2*Neq+4) = yb(4,2) - ya(4,3); % flux continuity between MPL & CL
% res(2*Neq+4) = yb(4,2); % zero flux at boundary
% res(2*Neq+4) = yb(3,2) - P_sat_o(T_L)/P; % saturated water vapor at liquid channel
res(4*Neq+3) = yb(3,2) - ya(3,3); % potential continuity between MPL & CL
res(4*Neq+4) = yb(4,3); % zero flux at boundary

% CARBON DIOXIDE GAS
res(0*Neq+5) = ya(5,1) - x_CO2_GC; % gas channel carbon dioxide content
res(0*Neq+6) = yb(6,1) - ya(6,2); % flux continuity between GDL & MPL
res(2*Neq+5) = yb(5,1) - ya(5,2); % potential continuity between GDL & MPL
res(2*Neq+6) = yb(6,2) - ya(6,3); % flux continuity between MPL & CL
% res(2*Neq+6) = yb(6,2) - i/(2*F); % zero flux at boundary
res(4*Neq+5) = yb(5,2) - ya(5,3); % potential continuity between MPL & CL
res(4*Neq+6) = yb(6,3); % zero flux at boundary

% CARBON MONOXIDE GAS
res(0*Neq+7) = ya(7,1) - x_CO_GC; % gas channel carbon monoxide content
res(0*Neq+8) = yb(8,1) - ya(8,2); % flux continuity between GDL & MPL
res(2*Neq+7) = yb(7,1) - ya(7,2); % potential continuity between GDL & MPL
res(2*Neq+8) = yb(8,2) - ya(8,3); % flux continuity between MPL & CL
% res(2*Neq+8) = yb(8,2) + i/(2*F); % zero flux at boundary
res(4*Neq+7) = yb(7,2) - ya(7,3); % potential continuity between MPL & CL
res(4*Neq+8) = yb(8,3); % zero flux at boundary

% TEMPERATURE
res(0*Neq+9) = ya(9,1) - T_C; % cathode boundary temperature
res(0*Neq+10) = yb(9,3) - T_L; % liquid boundary temperature
for d = 2:Nd
    res(2*(d-1)*Neq+ 9) = ya(9,d) - yb(9,d-1); % potential continuity
    res(2*(d-1)*Neq+10) = ya(10,d) - yb(10,d-1); % flux continuity
end
    
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

end