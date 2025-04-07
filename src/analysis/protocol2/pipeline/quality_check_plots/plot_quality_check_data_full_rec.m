function plot_quality_check_data_full_rec(f_data, v_data)

figure; 

subplot(5,1,1)
plot(f_data)
ax = gca;
ax.XAxis.Visible = 'off';
title('f-data')
xlim([0 numel(f_data)])

subplot(5,1,2:3)
plot(v_data)
xlim([0 numel(v_data)])
ylabel('v-data')
ax = gca;
ax.XAxis.Visible = 'off';
title(strcat("var - v-data = ", string(var(v_data))))

subplot(5,1,4:5)
plot(movmean(v_data, 20000))
xlim([0 numel(v_data)])
ylabel('movmean(v-data, 20000)')
hold on 
plot(movmean(v_data, 100000))

title(strcat("var - movmean - 100000 = ", string(var(movmean(v_data, 100000)))))
f = gcf;
f.Position = [25 629 1121 392];

end 













