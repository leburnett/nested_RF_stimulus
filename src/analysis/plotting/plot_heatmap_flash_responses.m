function f = plot_heatmap_flash_responses(data_comb)
% PLOT_HEATMAP_FLASH_RESPONSES  Create spatial heatmap of receptive field.
%
%   F = PLOT_HEATMAP_FLASH_RESPONSES(DATA_COMB) generates a 2D heatmap
%   showing the spatial distribution of neural responses across the
%   flash stimulus grid.
%
%   INPUT:
%     data_comb - NxN matrix of response values at each grid position
%                 N=14 for 4px flashes, N=10 for 6px flashes
%                 Positive values indicate excitation, negative = inhibition
%
%   OUTPUT:
%     f - Figure handle
%
%   FIGURE:
%     - Red/blue diverging colormap (redblue)
%     - Color limits centered on median value
%     - Square aspect ratio
%     - Colorbar showing response magnitude
%
%   PURPOSE:
%     Quick visualization of receptive field structure showing
%     spatial location of excitatory center and inhibitory surround.
%
%   See also PARSE_FLASH_DATA, GAUSSIAN_RF_ESTIMATE, REDBLUE

    figure; 
    imagesc(data_comb)
    cmap = redblue();
    colormap(cmap)
    
    med_val = median(data_comb(:));
    clim([med_val-0.5 med_val+0.5])
    colorbar
    axis square
    f = gcf;
end 