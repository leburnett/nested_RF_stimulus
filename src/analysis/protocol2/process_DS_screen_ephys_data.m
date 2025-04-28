function process_DS_screen_ephys_data()

exp_folder = cd;

PROJECT_ROOT = "/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2";

strain_str = "TmY3";

resultant_angle = process_bars_p2(exp_folder, strain_str, PROJECT_ROOT);

process_flash_p2(exp_folder, strain_str, PROJECT_ROOT, resultant_angle)


end 