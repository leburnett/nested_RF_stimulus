function generate_boxplot_data_tbl(data_tbl, data_type)

    % 1 - Add group column ('gp') to the results table:
    
    for i = 1:height(data_tbl)
    
        if data_tbl.strain(i) == "control"
            if data_tbl.type(i) == "OFF"
                data_tbl.gp(i) = 1; % control OFF
            else
                data_tbl.gp(i) = 2; % control ON
            end 
        else
            if data_tbl.type(i) == "OFF"
                data_tbl.gp(i) = 3; % ttl OFF
            else
                data_tbl.gp(i) = 4; % ttl ON
            end 
        end 
    end 
    
    % 2- Plot figure with boxchart of data
    
    figure
    b = boxchart(data_tbl.gp, data_tbl.(data_type), "GroupByColor", data_tbl.gp, "BoxWidth", 2.5);
    % Update colour of bars:
    b(2).BoxFaceColor = [0.3010 0.7450 0.9330];
    b(3).BoxFaceColor = [0.6050 0.0780 0.1340];
    b(4).BoxFaceColor = [0.9350 0.4780 0.6040];

    xticks([1,2,3,4])
    xticklabels({'ctrl-off', 'ctrl-on', 'ttl-off', 'ttl-on'})
    % xlabel('group')
    ylabel(strrep(data_type, '_', '-'))
    xlim([0 5])

    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.LineWidth = 1.2;
    ax.FontSize = 14;
    
    f = gcf;
    f.Position = [798   633   228   401];

end 