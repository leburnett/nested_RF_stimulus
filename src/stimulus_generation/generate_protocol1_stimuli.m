function generate_protocol1_stimuli(hemi)
% Generate both the pattern and position function for RF flash stimuli - protocol 1. 
% Non-overlapping grid of flashes that cover ~1/2 of the arena.
% Initial protocol to find the centre of the receptive field.
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
            params.root_dir = fullfile(ROOT_DIR, 'protocols', 'bkg4', 'RHS');
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
    t_bkg_s = 2; % time in seconds to display the background frame before beginning the flashes.
    generate_static_pos_fn(func_save_dir, t_bkg_s)
    
end 
