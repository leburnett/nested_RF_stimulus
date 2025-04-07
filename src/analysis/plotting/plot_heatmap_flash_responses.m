function f = plot_heatmap_flash_responses(data_comb)
    
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