function generate_protocol1_stimuli(hemi)
% GENERATE_PROTOCOL1_STIMULI  Generate Protocol 1 flash stimuli for initial RF mapping.
%
%   GENERATE_PROTOCOL1_STIMULI(HEMI) creates pattern and position function
%   files for Protocol 1, which presents a non-overlapping grid of flashes
%   covering half the arena to identify the approximate receptive field
%   location of the recorded neuron.
%
%   INPUT:
%     hemi - String specifying which hemisphere to stimulate:
%            'left'  - Presents stimuli on the fly's left visual field
%                      (columns 17-112, for recording from left hemisphere)
%            'right' - Presents stimuli on the fly's right visual field
%                      (columns 81-180, for recording from right hemisphere)
%
%   PROTOCOL 1 PURPOSE:
%     This is the first phase of the two-protocol RF mapping system.
%     Protocol 1 presents coarse flash grids (12px and 6px) to rapidly
%     identify which region of the screen elicits neural responses.
%     The experimenter observes the responses and identifies the "peak_frame"
%     (frame number showing the best response), which becomes the centerpoint
%     for the high-resolution Protocol 2.
%
%   STIMULUS PARAMETERS:
%     Flash sizes: 12px and 6px (both generated)
%     Overlap: 0% (non-overlapping grid)
%     Flash duration: 200ms
%     Inter-flash interval: 150ms
%     ON/OFF: Both contrasts presented
%     Pixel intensity: [bkg=4, off=0, on=15]
%     Background display: 5 seconds before flashes begin
%
%   OUTPUT:
%     Creates files in nested_RF_stimulus/protocols/bkg4/[LHS|RHS2]/:
%       - Patterns/ folder with .pat and .mat files
%       - Functions/ folder with .pfn and .mat files
%       - Log Files/ folder (created when experiment runs)
%
%   EXAMPLE:
%     % For recording from left hemisphere (stimulate right visual field)
%     generate_protocol1_stimuli('left')
%
%     % For recording from right hemisphere (stimulate left visual field)
%     generate_protocol1_stimuli('right')
%
%   See also GENERATE_STIMULUS, GENERATE_FLASH_PATTERN, GENERATE_FLASH_FUNCTION,
%            GENERATE_STATIC_POS_FN, GENERATE_PROTOCOL2

ROOT_DIR = 'C:\matlabroot\GitHub\nested_RF_stimulus';

    for flash_sz = [12, 6]
    
        %% Flash size in pixels
        params.flash_sz_px = flash_sz; % Size of flash in pixels
        
        %% Region of arena on which to present the stimulus. 
        assert(ismember(hemi, {'left', 'right'}), ...
        'Input argument "hemi" must be either "left" or "right".');
        
        if hemi == "left"
            %  1 - For the 1/2 screen area. 6 panels. RHS of screen when look from the
            %  front. LHS from fly's perspective. For recording from left hemisphere. 
            disp_x1 = 17;
            disp_x2 = 112;
            disp_y1 = 1;
            disp_y2 = 48;
            params.root_dir = fullfile(ROOT_DIR, 'protocols', 'bkg4', 'LHS');
        elseif hemi == "right"
            %  2 - For the 1/2 screen area. 6 panels. LHS of screen when look from the
            %  front. RHS from fly's perspective. For recording from right hemisphere. 
            % disp_x1 = 97;
            % disp_x2 = 196;
            % disp_y1 = 1;
            % disp_y2 = 48;
            % params.root_dir = fullfile(ROOT_DIR, 'protocols', 'RHS');
            disp_x1 = 81;
            disp_x2 = 180;
            disp_y1 = 1;
            disp_y2 = 48;
            params.root_dir = fullfile(ROOT_DIR, 'protocols', 'bkg4', 'RHS2');
        end 
        
        %% Other parameters which don't change with the area of the screen that the
        % pattern is being presented. 
       
        params.protocol = 'protocol1';
        params.px_intensity = [4, 0, 15]; % Intensity parameters of flash stimulus - [bkg_color, off_color, on_color]
        params.px_rng = [disp_y1, disp_y2, disp_x1, disp_x2]; % Pixel range of screen to present on - [row_start, row_end, col_start, col_end]
        params.overlap = 0; % Overlap of flashes - between 0 and 1. 
        
        params.bkg_frame = 1; % The background frame = frame 1 in the pattern:
        params.interval_dur = 0.15; % duration of interval background screen in seconds. Used to be 0.05.
        params.flash_dur = 0.2; % duration of flash in seconds.
        params.on_off = "both"; % For protocol 1, show both ON and OFF flashes.
        
        generate_stimulus(params)
    
    end 

    % Generate static position function for presenting background frame at
    % the beginning of the protocol:
    func_save_dir = fullfile(params.root_dir, 'Functions');
    t_bkg_s = 5; % time in seconds to display the background frame before beginning the flashes.
    generate_static_pos_fn(func_save_dir, t_bkg_s)
    
end 
