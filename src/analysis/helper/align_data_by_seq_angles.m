function data_ordered = align_data_by_seq_angles(data)
    % Rearrange the timeseries in 'data' so that they are arranged in
    % sequential angular order. i.e. 0, 1/16pi, 2/16pi etc.. 

    % Initially they are ordered in PD/ND pairs - so 0, pi, etc..

    % Creates 'data2' which is the same size as 'data'. Move all data for
    % individual reps as well as the mean.
    
    plot_order= [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

    data_ordered = cell(32, 4);

    n_dir_and_sp = height(data);
    n_speeds = 2;
    n_dir = n_dir_and_sp/n_speeds;

    for k = 1:n_dir_and_sp
   
        if k >n_dir % 56 dps
            data_idx = plot_order(k-n_dir);
            data_idx = data_idx+n_dir;
        else % 28 dps
            data_idx = plot_order(k);
        end 
        data_ordered(k, :) = data(data_idx, :);
    end 

end 




