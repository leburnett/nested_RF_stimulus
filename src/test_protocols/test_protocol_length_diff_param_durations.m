%% Calculate how long protocol 1 will last with different parameter values:

% These parameters are fixed:
n_12px_flashes = 32;
n_6px_flashes = 128;
n_on_off = 2;

% Parameters that can be changed:
t_interval_start = 2;
n_reps = 4;
t_flash = 0.2;
t_bkg = 0.1; % originally 0.05; 50ms.

% Calculate
total_t_secs = n_reps * (t_interval_start + n_on_off*((t_flash+t_bkg)*n_12px_flashes)+ n_on_off*((t_flash+t_bkg)*n_6px_flashes));
total_t_mins = total_t_secs/60;

% 5.4667 minutes for P1 (all reps) - 0.2 + 0.05 - 2s bkg before flashes
% 8.5667 minutes for P1 (all reps) - 0.2 + 0.2 
% 7.5 minutes for P1 (all reps) - 0.2 + 0.15 


%% Protocol 2: 

n_reps = 3; 
n_flashes = 196;

sflash_dur = 0.16;
sflash_bkg = 0.44; % was 0.34
slow_flash_t = sflash_dur + sflash_bkg;
total_sflash = slow_flash_t * n_flashes; % seconds

fflash_dur = 0.08;
fflash_bkg = 0.22; % was 0.17
fast_flash_t = fflash_dur + fflash_bkg;
total_fflash = fast_flash_t * n_flashes; % seconds

% Bar patterns are made using the GUI, then a function "crops" a 30px square around the
% centre of the screen and then position that square centred on the square
% found to give the largest response in P1.

% TODO - update to make these initial patterns from script too.... 
% Or, make the cropped and repositioned bar patterns entirely from
% script....?

% 1 - Increase the gap between moving bars - change the position function. 
% 2 - Make the bars thinner??? 

dur_sbar = 2.328;
dur_fbar = 1.164;
static_t = 1; % was 0
n_dir = 16; % 8 orientations - both directions. 

total_sbars = (dur_sbar + static_t) * n_dir; 
total_fbars = (dur_fbar + static_t) * n_dir; 

total_t_secs = (total_sflash + total_fflash + total_sbars + total_fbars) * n_reps; % [slow flashes, fast flashes, slow bars, fast bars]
total_t_mins = total_t_secs/60;

% 3.3812 minutes per rep. 
% 10.1436 minutes for all 3 reps. 

% With 1s static before each bar direction
% 11.7436 minutes

% Check the background pixel value of the bars.... seems brighter than the
% normal background... 

