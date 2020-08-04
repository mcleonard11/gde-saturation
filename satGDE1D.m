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
RH = 0.90; % [-] relative humidity in gas channel
s_C = 0.01; % [-] liquid water saturation at GDL/GC interface
T = T_0+30; % [K] temperature of gas channel

N_w = [0.5 1 1.5 2.0]'; % water injection flux (proxy for water generation)

% MATERIAL PARAMETERS
L = [100 10]*1e-6; % [m] gas diffusion layer thicccness
s_im = s_C; % [-] immobile liquid water saturation
eps_p_GDL = 0.76; % [-] porosity of GDL
eps_p_CL = 0.4; % [-] porosity of CL
kappa_GDL = 6.15e-12; % [m^2] absolute permeability of GDL
kappa_CL = 1e-13; % [m^2] absolute permeability of CL
tau_GDL = 1.6; % [-] pore tortuosity of GDL
tau_CL = 1.6; % [-] pore tortuosity of CL

% WATER CONSTITUTIVE RELATIONSHIPS
P_sat = @(T) exp(23.1963-3816.44./(T-46.13)); % [Pa] saturation pressure of water vapor
mu = @(T) 1e-3*exp(-3.63148+542.05./(T-144.15)); % [Pa*s] dynamic viscosity of liquid water

% MODEL PARAMETERIZATION
D = @(eps_p,tau,s,P,T) eps_p/tau^2*(1-s).^3.*(T/T_ref).^1.5*(P_ref/P); % [-] scaling factor for gas diffusivities
D_H2O = @(eps_p,tau,s,T) 0.36e-4*D(eps_p,tau,s,P,T); % [m^2/s] H2O gas phase diffusion coefficient
x_H2O_C = RH*P_sat(T)/P; % [-] mole fraction of water vapor in gas channel

% AUXILIARY FUNCTIONS
iff = @(cond,a,b) cond.*a + ~cond.*b; % vectorized ternary operator

% MATERIAL CONSTITUTIVE RELATIONSHIPS
s_red = @(s) (s-s_im)/(1-s_im); % reduced liquid water saturation
gamma_ec = @(x_H2O,x_sat,s,T) 2e6*iff(x_H2O<x_sat,5e-4*s_red(s),6e-3*(1-s_red(s))).*sqrt(R*T/(2*pi*M_w)); % [1/s] evaporation/condensation rate
dpds = @(s) 0.00011*44.02*exp(-44.02*(s-0.496))+278.3*8.103*exp(8.103*(s-0.496)); % [Pa] derivative of capillary pressure-saturation relationship of GDL
D_s = @(kappa,s,T) kappa*(1e-6+s_red(s).^3)./mu(T).*dpds(s); % [m^2/s] liquid water transport coefficient


% INITIAL MESH
Lsum = [0 cumsum(L)];
Nd = numel(L); % number of domains
x = interp1(0:Nd, Lsum, linspace(0, Nd, Nd*10+1));
x = sort([x Lsum(2:end-1)]); % duplicate interface nodes

% SOLVER PREPARATION
sol = bvpinit(x, @yinit);
options = bvpset('Vectorized', 'on', 'NMax', 1e3, 'RelTol', 1e-4, 'AbsTol', 1e-6);

% PARAMETER SWEEP
SOL = cell(size(N_w));
Np = numel(N_w); % number of parameters in the sweep
Neq = size(sol.y,1)/2; % number of 2nd-order differential equations
for k = 1:Np
    sol = bvp4c(@odefun, @(ya,yb) bcfun(ya, yb, N_w(k)), sol, options);
    SOL{k} = sol;
end
% POSTPROCESSING
Nref = 2; % number of refinements for smoother curve plotting
domains = [1 1;
           1 1];
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
unit_scale = [1 1;
              1 1];
quantity = {'s','x_{H_2O}';
            'j_s','j_{H2O}'};
c = winter(Np);
for m = 1:2
    figure('Name', fig_names{m})
    for n = 1:Neq
        subplot(1,2,n)
        box on
        hold on
        us = unit_scale(m,n);
        for k = 1:Np
            plot(SOL{k}.x*1e6, SOL{k}.y(2*(n-1)+m,:)*us, 'Color', c(k,:), 'DisplayName', [num2str(N_w(k)) ' blah'])
        end
        xlim([Lsum(find(domains(n,:),1,'first')) Lsum(find(domains(n,:),1,'last')+1)]*1e6)
        ylim(ylim)
        xlabel('x [um]')
        ylabel(quantity(m,n))
        for x = Lsum(2:end-1)
            l = line([x x]*1e6, ylim, 'Color', 'k');
            set(get(get(l, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
        end
    end
    legend(cellstr(num2str(N_w)),'Location','best');
end

function dydx = odefun(x, y, subdomain)

% READ POTENTIALS & FLUXES
s   = y( 1,:); j_s   = y( 2,:);
x_w = y( 3,:); j_x_w = y( 4,:);

% ZERO-INITIALIZE ALL DERIVATIVES
z = zeros(size(x));
ds    = z; dj_s    = z;
dx_w = z; dj_x_w = z;

% COMPUTE DERIVATIVES
switch subdomain
    case 1 % GAS DIFFUSION LAYER
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T)./P; % saturation water vapor mole fraction
        S_ec = gamma_ec(x_w,x_sat,s,T).*C.*(x_w-x_sat); % evaporation/condensation reaction rate
        ds = -j_s./((rho_w/M_w)*D_s(kappa_GDL,s,T)); % liquid water flux: j_s = -rho_w*D_s*grad(s) 
        dx_w = -j_x_w./(C.*D_H2O(eps_p_GDL,tau_GDL,s,T)); % water vapor flux: j_H2O_v = -rho_g*D_H2O*grad(x_H2O)
        dj_s = S_ec; % conservation of liquid water: div(j_s) = S_s
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_v) = S_H2O
    case 2 % CATALYST LAYER
        C = P./(R*T); % gas phase density
        x_sat = P_sat(T)./P; % saturation water vapor mole fraction
        S_ec = gamma_ec(x_w,x_sat,s,T).*C.*(x_w-x_sat); % evaporation/condensation reaction rate
        ds = -j_s./((rho_w/M_w)*D_s(kappa_CL,s,T)); % liquid water flux: j_s = -rho_w*D_s*grad(s)
        dx_w = -j_x_w./(C.*D_H2O(eps_p_CL,tau_CL,s,T)); % water vapor flux: j_H2O_v = -rho_g*D_H2O*grad(x_H2O)
        dj_s = S_ec; % conservation of liquid water: div(j_s) = S_s
        dj_x_w = -S_ec; % conservation of water vapor: div(j_H2O_v) = S_H2O
end

% ASSEMBLE DERIVATIVES
dydx = [ds   ; dj_s  ;
        dx_w ; dj_x_w];
end

function y0 = yinit(x, subdomain)
% POTENTIALS INITIALLY SET TO GAS CHANNEL CONDITIONS
s = s_C;
x_w = x_H2O_C;

% ALL FLUXES ARE INITIALLY ZERO
y0 = [s   ; 0; 
      x_w ; 0];
end

function res = bcfun(ya, yb, N_w)

res = ya(:); % homogeneous BC everywhere by default

% LIQUID WATER
res(0*Neq+1) = ya(1,1) - s_C; % GC liquid water content
res(0*Neq+2) = yb(2,1) - ya(2,2); % flux continuity between GDL & CL  
res(2*Neq+1) = ya(1,2) - yb(1,1); % potential continuity between GDL & CL
res(2*Neq+2) = yb(2,2) + N_w; % flux at CL water interface

% WATER VAPOR
res(0*Neq+3) = ya(3,1) - x_H2O_C; % gas channel water vapor content
res(0*Neq+4) = yb(4,1) - ya(4,2); % flux continuity between GDL & CL
res(2*Neq+3) = yb(3,1) - ya(3,2); % potential continuity between GDL & CL
res(2*Neq+4) = yb(4,2); % zero flux at boundary
    
end


end