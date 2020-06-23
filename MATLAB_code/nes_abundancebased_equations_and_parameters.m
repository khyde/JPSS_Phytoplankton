%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NES abundance-based model equations and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------
% NES Brewin model (B-NES)
% ------------------------

% Parameters:
% -----------
Cmpn = 0.8147; 
Dpn = 0.7776;
Cmp = 0.1509; 
Dp = 0.5429; 

% Equations (only input is satellite chl, denoted C_sat [not log-transformed!]):
% ------------------------------------------------------------------------------

% Size fractions (0-1)
Fpico_bnes = (Cmp .* (1 - exp(-1 .* (Dp ./ Cmp) .* C_sat))) ./ C_sat;
Fpiconano_bnes = (Cmpn .* (1 - exp(-1 .* (Dpn ./ Cmpn) .* C_sat))) ./ C_sat; %"piconano" refers to the combined pico and nano populations
Fnano_bnes = Fpiconano_bnes - Fpico_bnes;
Fmicro_bnes = (C_sat - (Cmpn .* (1 - exp(-1 .* (Dpn ./ Cmpn) .* C_sat)))) ./ C_sat;

% Size-specific chl [mg/m^3]
Cpico_bnes = Fpico_bnes .* C_sat;
Cpiconano_bnes = Fpiconano_bnes .* C_sat;
Cnano_bnes = Fnano_bnes .* C_sat;
Cmicro_bnes = Fmicro_bnes .* C_sat;

% --------------------------------
% NES SST Brewin model (B-NES-SST)
% --------------------------------

% Parameters: Same as above, but found based on closest matching SST in the
% LUT (contained in "nes_brewin_look_up_table_SST.mat")

% Equations: Same as above, but using the SST-dependent parameters.

% ------------------------
% NES Hirata model (H-NES)
% ------------------------

% Parameters:
% -----------
b1_micro = 1.0297;
b2_micro = -1.6841;
b3_micro = -0.1216;
b1_pico = -3.4459;
b2_pico = 0.6681;
b3_pico = 2.2859;

% Equations (input for size fractions is log10(C_sat), denoted logC_sat):
% -----------------------------------------------------------------------

% Size fractions (0-1)
Fpico_hnes = 1 ./ ((b1_pico) + exp(b2_pico .* logC_sat + b3_pico));
Fpico_hnes(Fpico_hnes>1) = 1; % correct for if fraction exceeds 1 or is less than 0
Fpico_hnes(Fpico_hnes<0) = 0;  
Fmicro_hnes = 1 ./ ((b1_micro) + exp(b2_micro .* logC_sat + b3_micro));
Fmicro_hnes(Fmicro_hnes>1) = 1;
Fmicro_hnes(Fmicro_hnes<0) = 0;  
Fnano_hnes = 1 - Fmicro_hnes - Fpico_hnes;
Fpiconano_hnes = Fpico_hnes + Fnano_hnes;

% Size-specific chl [mg/m^3]
Cmicro_hnes = Fmicro_hnes .* C_sat; % C_sat not log-transformed
Cnano_hnes = Fnano_hnes .* C_sat;
Cpico_hnes = Fpico_hnes .* C_sat;
Cpiconano_hnes = Fpiconano_hnes .* C_sat;
