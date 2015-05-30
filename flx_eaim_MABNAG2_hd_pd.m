function flx = flx_eaim_MABNAG2_hd_pd(time, input, R,k,T,press, M, rhoi,...
    p_sat, sigmai, D, alpham, NumbConc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Function for calculating mass fluxes of acid to the particle.           %
% This version is for system of 2 acids, 2 bases and water.               %
% 
% flx --> Particle Phase Fluxes [kg m-3 s-1] (positive is growth)
%           1) Water
%           2) SO4--
%           3) Organic Acid
%           4) NH3
%           5) Amine
%         Gas Phase Fluxes [kg m-3 s-1] (positive is growth)
%           6) Water
%           7) H2SO4
%           8) DMA
%           9) NH3
%          10) Amine
%                                                                         %
% This version (_pd) is for calculating the mass flux according to        %
% equations where particle diffusion and vapor molecule dimensions are    %
% taken into account according to Lehtinen & Kulmala 2003 (ACP) and       %
% Nieminen et al. 2010 (ACP).                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Avogadro constant (1/mol)
NA = 6.0221415e+23;  

% Use Masses in Gas Phase to Calculate Partial Pressures
pvi = input(6:10) ./ transpose( M(6:10) ) .* (R.*T);  %[Pa]

% Masses of each component [water acid1 acid2 base1 base2]
mpi = transpose( input(1:5) ) ./ NumbConc;  %[kg/particle]

% Particle mass [kg]
mp = sum(mpi);

% Mass fractions
Xmass = mpi./mp;

% Mole fractions
Xmole = (Xmass ./ M(1:5) )./sum(Xmass./ M(1:5) );

% Particle density [kg m-3]
rho = ( sum(Xmass./rhoi) )^(-1);

% Particle radius [m]
rp = (3*(mp/rho)/(4*pi)).^(1.0/3.0);

% Particle surface tension
sigma = sum(Xmole.*sigmai);

% Kelvin effect for the acids
Ke = exp(2.* M(1:5) .*sigma./R./T./rp./rho);

% Calculate Equilibrium Vapor Pressure as Saturation Pressure times Kelvin
% Term
p_eq = p_sat .* Ke;

% Particle diameter
dp = 2.0*rp;


%--- Transition regime corrections ----------------------------------------
% According to Lehtinen & Kulmala and Nieminen et al.

% Calculating diffusion coefficient of particle
visc_air = 1.84e-5; %Gas viscosity of air (kg/m/s)
M_air = 29.0e-3;  %Molar mass of air (kg/mol)
lambda_air = 2*visc_air/(press*(8*M_air/(pi*R*T))^(1/2));  %Mean free bath of air molecules
Kn_air = 2*lambda_air/dp;  %Knudsen number using air mean free bath
Cc = 1+Kn_air*(1.257+0.4*exp(-1.1/Kn_air)); %Cunningham slip correction factor
Diffp = k*T*Cc/(3*pi*visc_air*dp);  %Diffusion coefficient of particle

% Mass of vapor molecules [kg molec-1]
mv = M(6:10)./NA;

% Mean thermal speed [m s-1]
cv = (8*k*T./(pi.*mv)).^(1/2);  %vapors
cp = (8*k*T/(pi*mp))^(1/2);     %particle

% Mean free path  [m]
lambda = 3.0.*(Diffp+D)./(cp.^2+cv.^2).^(1/2);  

%Diameter of vapor molecules [m]
dv = (6.*mv./(pi.*rhoi)).^(1/3);

%Knudsen number
Kn = 2.*lambda./(dp+dv); 

%Transition regime correction factor for mass flux
beta = (1.0 + Kn)./(1.0 + (4./(3.0.*alpham) + 0.377).*Kn + 4.0./(3.0.*alpham).*Kn.^2.);

%--------------------------------------------------------------------------    


%--- Mass fluxes ----------------------------------------------------------
for j = 1:5
    if p_eq(j) < 1e-15
        p_eq(j) = 0.;
    end
    
    if (mpi(j) <= 0.0 & (pvi(j)-p_eq(j)) < 0)
        flx_mass(j) = 0.0;
    else
        flx_mass(j) = NumbConc .* M(j) * 2*pi*(dv(j) + dp)*(D(j) + Diffp) ...
                      * beta(j)*(pvi(j)-p_eq(j))./(R*T);
    end
end


% Pass Mass Fluxes to Function Output
flx(1:5,1)  =  flx_mass;    %To/from Particle Phase [kg m-3 s-1]
flx(6:10,1) = -flx_mass;  %To/from Vapor Phase [kg m-3 s-1] 

%--------------------------------------------------------------------------

end