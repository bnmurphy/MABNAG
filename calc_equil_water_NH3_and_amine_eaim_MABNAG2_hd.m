% changed by LA
%function [p_eq,moles_water,moles_NH3,moles_amine,radius,moles_ions] = calc_equil_water_NH3_and_amine_eaim_MABNAG2_hd(input_moles,RH,pvi,T,M,rhoi,sigmai,NumbConc)
function [p_sat] = calc_equil_water_NH3_and_amine_eaim_MABNAG2_hd(...
        input_moles,RH,T,M,rho,sigma,NumbConc,rp)
% end LA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Function for calculating equilibrium water, NH3 and amine content of    %
% a particle using E-AIM.                                                 %
%                                                                         %
% This is for system with 5 compounds:                                    %
% Compound 1 = water                                                      %
% Compound 2 = acid 1 (sulfuric acid)                                     %
% Compound 3 = acid 2 (organic acid)                                      %
% Compound 4 = base 1 (ammonia)                                           %
% Compound 5 = base 2 (amine)                                             %
%                                                                         %
% Input: input_moles = moles of each of the compound in a particle        %
%        RH = relative humidity                                           %
%        pvi = partial pressures of the compounds in the gas phase        %
%              [H2O H2SO4 org.acid NH3 amine]                             %
%        T = temperature                                                  %
%        M = molar masses [H2O H2SO4 org.acid NH3 amine]                  %
%        rhoi = liquid desities [H2O H2SO4 org.acid NH3 amine]            %
%        sigmai = surface tension [H2O H2SO4 org.acid NH3 amine]          %
%        NumbConc = particle number concentration (particles/m^-3)        %
%                                                                         %
% Result: p_eq = equilibrium vapor pressures of the compounds. Order of   %
%                compounds:                                               %
%                1 = water                                                %
%                2 = acid 1 (sulfuric acid)                               %
%                3 = acid 2 (organic di-acid)                             %
%                4 = base 1 (ammonia)                                     %
%                5 = base 2 (amine)                                       %
%                6 = neutral organic acid (NaN here as Neutral organic    %
%                    acid is not included)                                %
%         moles_water = number of moles of water in the particle in       %
%                       equilibrium                                       %
%         moles_NH3 = number of moles of NH3 in the particle in           %
%                     equilibrium                                         %
%         moles_amine = number of moles of amine in the particle in       %
%                       equilibrium                                       %
%         radius = radius of the particle                                 %
%         moles_ions = moles of each the ions and molecules in the        %
%                      particle. Order of the compounds:                  %
%                      1 = nH(aq)                                         %
%                      2 = nNH4(aq)                                       %
%                      3 = nDMA+(aq)                                      %
%                      4 = nHSO4(aq)                                      %
%                      5 = nSO4(aq)                                       %
%                      6 = nOH(aq)                                        %
%                      7 = nMalo-(aq)                                     %
%                      8 = nMalo2-(aq)                                    %
%                      9 = nH2O(aq)                                       %
%                      10 = nNH3(aq)                                      %
%                      11 = nMaloni(aq)                                   %
%                      12 = nDMA(aq)                                      %
%                      13 = nNeut_o(aq) (NaN here as Neutral organic is   %
%                           not included)                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Gas constant (J K^-1 mol^-1)
R = 8.314472;
% Boltzmann constant (m^2 kg s^-2 K^-1)
kb = 1.3806503e-23;
% Avogadro constant (mol^-1)
NA = 6.0221415e+23;


% Moles per m^3 in particle phase and vapor phase
moles = input_moles;  %[mol/m3]

radius = rp; %

%Kelvin effect based on mass weighted density and mole weighted surface
%tension
Ke = exp(2.*M.*sigma./R./T./radius./rho);

%Moles of NH3 per m^3 in gas phase 
n_g_NH3 = moles(9);

%Moles of amine per m^3 in gas phase 
n_g_amine = moles(10);

%Moles of Water per m^3 in gas phase
n_g_H2O = moles(6);
%Total moles of water per m^3
n_t_H2O = n_g_H2O + moles(1);  %[mol/m3]

%**********************************************************************
% Generating input file for E-AIM                                     *
%**********************************************************************

input_eaim = [298.15  1.0  1.0   1  1  0.42 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0   3  3  0  4  3 17  1 4  3 4  4 4  5 4  6 4  7 4  8 4  9 4  10 4  11 4  12 4  13 4  14 4  15 4  16 4  28 4  29 4  0.0  0.0  4  4  3  3  3  3  0  0  3  3];

% Columns:
% 1 = Temperature (K)
input_eaim(1,1) = T;
% 2 = system pressure (atm)
input_eaim(1,2) = 1.0;
% 3 = system volume (m3)
input_eaim(1,3) = 1.0;
% 4 = 1 = relative humidity of the system is to be fixed, 2  = total number of moles of water in the system is given., 3 = total number of moles of water in the system, but water is not allowed to partition into the vapour
input_eaim(1,4) = 2;  %Used to be 1
% 5 = Water dissociation optio: -1 = dissociation not calculated, 0 = dissociation calculated if either H+ of OH- is given as input, 1 = water dissociation calculated always
input_eaim(1,5) = 1;
% 6 = RH (as fraction (0.1 < RH <0.99)) or total number of moles of water in the system
input_eaim(1,6) = n_t_H2O;    %moles of water in the system; used to be RH
% 7 = Number of moles of H+
% 8 = Number of moles of NH4+
% 9 = Number of moles of Na+ (must be zero for model II)
% 10 = Number of moles of SO42-
% 11 = Number of moles of NO3-
% 12 = Number of moles of Cl- (must be zero for model II)
% 13 = Number of moles of Br- (must be zero for model II)
% 14 = Number of moles of OH-
% 15 = Number of moles of NH3
% 16 = control parameter for gas phase HNO3. 0 = partitioning to gas phase is allowed, 3 = partitioning to gas phase is NOT allowed, 4 = partitioning to gas phase is NOT allowed but the equilibrium vapor pressure over liquid is reported
input_eaim(1,16) = 3;
% 17 = control parameter for gas phase HCl. 0 = partitioning to gas phase is allowed, 3 = partitioning to gas phase is NOT allowed, 4 = partitioning to gas phase is NOT allowed but the equilibrium vapor pressure over liquid is reported (must be zero for model II)
input_eaim(1,17) = 3;
% 18 = control parameter for gas phase NH3. 0 = partitioning to gas phase is allowed, 3 = partitioning to gas phase is NOT allowed, 4 = partitioning to gas phase is NOT allowed but the equilibrium vapor pressure over liquid is reported
input_eaim(1,18) = 0;
% 19 = control parameter for gas phase H2SO4. 0 = partitioning to gas phase is allowed, 3 = partitioning to gas phase is NOT allowed, 4 = partitioning to gas phase is NOT allowed but the equilibrium vapor pressure over liquid is reported
input_eaim(1,19) = 0; %Use to be 4
% 20 = control parameter for gas phase HBr. 0 = partitioning to gas phase is allowed, 3 = partitioning to gas phase is NOT allowed, 4 = partitioning to gas phase is NOT allowed but the equilibrium vapor pressure over liquid is reported (must be zero for model II)
input_eaim(1,20) = 3;
% 21 = number of solids whose options are to be individually entered. NOTE! here zero always, chack why e-aim does not accept the solid options
input_eaim(1,21) = 17;
% 22-55 solid options
% 56 = the number of moles of the organic compound 1. (org. acid)
% 57 = the number of moles of the organic compound 2. (amine)
% 58 = Gas option for organic species 1. 0 = partitioning to gas phase is allowed, 3 = partitioning to gas phase is NOT allowed, 4 = partitioning to gas phase is NOT allowed but the equilibrium vapor pressure over liquid is reported.
input_eaim(1,58) = 0;  %Use to be 4
% 59 = Gas option for organic species 2. 0 = partitioning to gas phase is allowed, 3 = partitioning to gas phase is NOT allowed, 4 = partitioning to gas phase is NOT allowed but the equilibrium vapor pressure over liquid is reported.
input_eaim(1,59) = 0;
% 60 = Solid option for organic species 1. 0 = a solid can form, 3 = the solid is not allowed to form, 4 = the solid is not be allowed to form and the degree of supersaturation will be reported
input_eaim(1,60) = 3;
% 61 = Solid option for organic species 2. 0 = a solid can form, 3 = the solid is not allowed to form, 4 = the solid is not be allowed to form and the degree of supersaturation will be reported
input_eaim(1,61) = 3;
% 62 = mixed solid options for organic species 1
input_eaim(1,62) = 3;
% 63 = mixed solid options for organic species 2
input_eaim(1,63) = 3;
% 64 = dissociation options for organic species 1
input_eaim(1,64) = 0;
% 65 = dissociation options for organic species 2
input_eaim(1,65) = 0;
% 66 = liquid/liquid equilibrium options for organic species 1
input_eaim(1,66) = 3;
% 67 = liquid/liquid equilibrium options for organic species 2
input_eaim(1,67) = 3;

%%%%%%%%% added by LA %%%%%%%%%%

input_eaim(1,68) = 0;
input_eaim(1,69) = 0;
%input_eaim
%%%%%%%%%% end LA %%%%%%%%%%%%%%%%%


% For H+, SO4 and org. acid, input the total moles:

% H+ and SO4:
if moles(2) < 0.1e-23
    input_eaim(1,7) = 2*0.1e-23;    %H+
    input_eaim(1,10) = 0.1e-23;     %SO4--
else
    input_eaim(1,7) = 2*( moles(2) + moles(7) );  %H+
    input_eaim(1,10) = moles(2) + moles(7);   %SO4--
end

% Organic acid:
if moles(3) < NumbConc*1.0e-26;
    input_eaim(1,56) = NumbConc*1.0e-26;
else
    input_eaim(1,56) = moles(3) + moles(8);
end

% For NH3 and amine number of moles per m^3 in gas phase is known.
% Input have to be the total number of moles per m^3 (=aqueous + gas phase)

%NH3:
n_g_NH3_eff = n_g_NH3;  %Effective gas phase concentration. This is used here since E-AIM calculates the equilibrium for a flat surface
input_eaim(1,15) = moles(4) + n_g_NH3_eff;   %total number of moles (= aqueous+liquid phase)

%Amine:
n_g_amine_eff = n_g_amine;  %Effective gas phase concentration. This is used here since E-AIM calculates the equilibrium for a flat surface
input_eaim(1,57) = moles(5) + n_g_amine_eff;   %total number of moles (= aqueous+liquid phase)


%Do File Structure Work
delete('eaim.dat')
copyfile('eaim_simon_LA.dat','eaim.dat')

%Write the E-AIM input values in a file
fid = fopen('eaim.dat','a');
%fprintf(fid, '%4.2f %4.1f %4.1f %1i %2i %6.4E %6.4E %2.2E %2.2E %6.4E %2.2E %2.2E %2.2E %2.2E %6.4E %1i %1i %1i %1i %2i %3i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %6.4E %6.4E %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i\n', input_eaim');
fprintf(fid, '%4.2f %4.1f %4.1f %1i %2i %6.4E %6.4E %2.2E %2.2E %6.4E %2.2E %2.2E %2.2E %2.2E %6.4E %1i %1i %1i %1i %2i %3i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %6.4E %6.4E %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i\n', input_eaim');
fclose(fid);

% Run E-AIM
system('eaimIntel.exe');


%**********************************************************************
% Read the results from E-AIM from the file saved by the .exe         *
%**********************************************************************

fid = fopen('eaim.r2a','r');
head = textscan(fid,'%s',50);
output_eaim = textscan(fid,'%f %f %c %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
%fid = fopen('test.res','r')
%output_eaim = fscanf(fid,' Output record  %i: T=%f  RH=%f nH(aq)= %E fH(aq)= %E nNH4(aq)= %E fNH4(aq)= %E nDMA+(aq)= %E fDMA+(aq)= %E nHSO4(aq)= %E fHSO4(aq)= %E nSO4(aq)= %E fSO4(aq)= %E nOH(aq)= %E fOH(aq)= %E nMalo-(aq)= %E fMalo-(aq)= %E nMalo2-(aq)= %E fMalo2-(aq)= %E nH2O(aq)= %E fH2O(aq)= %E nNH3(aq)= %E fNH3(aq)= %E nMaloni(aq)= %E fMaloni(aq)= %E nDMA(aq)= %E fDMA(aq)= %E pH2O(g)= %E pNH3(g)= %E pH2SO4(g)= %E pMaloni(g)= %E pDMA(g)= %E\n');
fclose(fid);

head = head{1};

head(3) = []; % remove Err header
if strmatch('R',output_eaim(3)) == 1    % 'R' (error) is sometimes included as en extra column on the 2nd row without changing number of columns in the header
    output_eaim(3) = [];
end

head(1:7) = [];   % removing columns to the left of T (character element must at least be removed in order to apply cell2mat
output_eaim(1:7) = [];
out_eaim = cell2mat(output_eaim);

if isempty(out_eaim) == 1
    fid = fopen('eaim.r2a','r');
    head = textscan(fid,'%s',50);
    output_eaim = textscan(fid,'%f %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
    %fid = fopen('test.res','r')
    %output_eaim = fscanf(fid,' Output record  %i: T=%f  RH=%f nH(aq)= %E fH(aq)= %E nNH4(aq)= %E fNH4(aq)= %E nDMA+(aq)= %E fDMA+(aq)= %E nHSO4(aq)= %E fHSO4(aq)= %E nSO4(aq)= %E fSO4(aq)= %E nOH(aq)= %E fOH(aq)= %E nMalo-(aq)= %E fMalo-(aq)= %E nMalo2-(aq)= %E fMalo2-(aq)= %E nH2O(aq)= %E fH2O(aq)= %E nNH3(aq)= %E fNH3(aq)= %E nMaloni(aq)= %E fMaloni(aq)= %E nDMA(aq)= %E fDMA(aq)= %E pH2O(g)= %E pNH3(g)= %E pH2SO4(g)= %E pMaloni(g)= %E pDMA(g)= %E\n');
    fclose(fid);
    
    head = head{1};
    
    head(3) = []; % remove Err header
    if strmatch('GR',char(output_eaim{3})) == 1    % 'R' (error) is sometimes included as en extra column on the 2nd row without changing number of columns in the header
        output_eaim(3) = [];
    end
    
    head(1:7) = [];   % removing columns to the left of T (character element must at least be removed in order to apply cell2mat
    output_eaim(1:7) = [];
    out_eaim = cell2mat(output_eaim);
    
end

%%%%%%% end LA %%%%%%%%%

%out_eaim = (output_eaim(1:32))';
% 1 = Output record
% 2 = T
% 3 = RH
% 4 = nH(aq)
% 5 = fH(aq)
% 6 = nNH4(aq)
% 7 = fNH4(aq)
% 8 = nDMA+(aq)
% 9 = fDMA+(aq)
% 10 = nHSO4(aq)
% 11 = fHSO4(aq)
% 12 = nSO4(aq)
% 13 = fSO4(aq)
% 14 = nOH(aq)
% 15 = fOH(aq)
% 16 = nMalo-(aq)
% 17 = fMalo-(aq)
% 18 = nMalo2-(aq)
% 19 = fMalo2-(aq)
% 20 = nH2O(aq)
% 21 = fH2O(aq)
% 22 = nNH3(aq)
% 23 = fNH3(aq)
% 24 = nMaloni(aq)
% 25 = fMaloni(aq)
% 26 = nDMA(aq)
% 27 = fDMA(aq)
% 28 = pH2O(g)
% 29 = pNH3(g)
% 30 = pH2SO4(g)
% 31 = pMaloni(g)
% 32 = pDMA(g)


moles_ions_out = [];


%%%%%%%% changed by LA %%%%%%

%moles_ions_out = [out_eaim(4) out_eaim(6) out_eaim(8) out_eaim(10) out_eaim(12) out_eaim(14) out_eaim(16) out_eaim(18) out_eaim(20) out_eaim(22) out_eaim(24) out_eaim(26) NaN];


moles_ions_out = [out_eaim(strmatch('n_H(aq)',head)) out_eaim(strmatch('n_NH4(aq)',head))...
    out_eaim(strmatch('n_DMA+(aq)',head)) out_eaim(strmatch('n_HSO4(aq)',head))...
    out_eaim(strmatch('n_SO4(aq)',head)) out_eaim(strmatch('n_OH(aq)',head))...
    out_eaim(strmatch('n_Malo-(aq)',head)) out_eaim(strmatch('n_Malo2-(aq)',head))...
    out_eaim(strmatch('n_H2O(aq)',head)) out_eaim(strmatch('n_NH3(aq)',head))...
    out_eaim(strmatch('n_Maloni(aq)',head)) out_eaim(strmatch('n_DMA(aq)',head)) NaN];

% 1 = nH(aq)
% 2 = nNH4(aq)
% 3 = nDMA+(aq)
% 4 = nHSO4(aq)
% 5 = nSO4(aq)
% 6 = nOH(aq)
% 7 = nMalo-(aq)
% 8 = nMalo2-(aq)
% 9 = nH2O(aq)
% 10 = nNH3(aq)
% 11 = nMaloni(aq)
% 12 = nDMA(aq)
% 13 = nNeut_o(aq) (NaN here as Neutral organic is not included)


%n_water = (out_eaim(14) + out_eaim(20));
%n_aq_NH3 = (out_eaim(6) + out_eaim(22));
%n_gas_NH3 = out_eaim(29)*101325/(R*T);
%n_aq_amine = (out_eaim(8) + out_eaim(26));
%n_gas_amine = out_eaim(32)*101325/(R*T);

fid = fopen('eaim.r3a','r');
output3a_eaim = textscan(fid,'%f',59,'headerlines',1,'delimiter','\t');
fclose(fid);
out3a_eaim = output3a_eaim{1}';

n_water = (out_eaim(strmatch('n_OH(aq)',head)) + out_eaim(strmatch('n_H2O(aq)',head)));
n_aq_NH3 = (out_eaim(strmatch('n_NH4(aq)',head)) + out_eaim(strmatch('n_NH3(aq)',head)));
n_gas_NH3 = out3a_eaim(55)*101325/(R*T);
n_aq_amine = (out_eaim(strmatch('n_DMA+(aq)',head)) + out_eaim(strmatch('n_DMA(aq)',head)));
n_gas_amine = out3a_eaim(59)*101325/(R*T);


%Equilibrium vapor pressures (Over flat surface)
%This is the main interest of Ben and Jan
p_sat = [out3a_eaim(1,52) out3a_eaim(1,56) out3a_eaim(1,58) ...
    out3a_eaim(1,55) out3a_eaim(1,59)].*101325;


moles_water = n_water/NumbConc;
moles_NH3 = n_aq_NH3/NumbConc;
moles_amine  = n_aq_amine/NumbConc;
moles_ions = moles_ions_out./NumbConc;
radius = rp;
moles_moles_NH3 = moles(4);
moles_moles_amine = moles(5);
moles_moles_water = moles(1);

