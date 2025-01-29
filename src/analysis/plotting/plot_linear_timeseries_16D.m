function plot_linear_timeseries_16D(data_ordered)

    figure;

    for i2 = 1:32

        subplot(2,16,i2)

        plot(data_ordered{i2, 1}, 'k');
        n_xvals = numel(data_ordered{i2, 1});
        hold on 
        plot([0 n_xvals], [median_voltage median_voltage], 'r', 'LineWidth', 0.3)

        box off
        ax = gca;
        % ax.XAxis.Visible = 'off';
        ax.TickLength = [0.04 0.04];
        ax.TickDir = 'out';

        if i2>16
            ylim([-75 -40])
            ax.XTickLabel = {'0', '0.5', '1'};
        else
            ylim([-75 -30])
            ax.XTickLabel = {'0', '1', '2'};
        end 

        if i2 ==1 || i2 == 17
            ylabel('Voltage (mV)');
            xlabel('Time (s)')
        else
            ax.YAxis.Visible = 'off';
        end 

        if i2 == 8 
            title('28 dps')
        elseif i2 == 24
            title('56 dps')
        end 

    end 

    f = gcf;
    f.Position = [5 743 1795 283];

end 