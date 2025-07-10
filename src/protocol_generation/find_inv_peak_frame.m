function find_inv_peak_frame(peak_frame)
% Find the 'peak_frame' value for the inverse contrast flash at the exact
% same position to 'peak_frame'. Only smaller, 6 pixel squares. 

if peak_frame < 130
    inv_index = peak_frame+128;
else
    inv_index = peak_frame-128;
end 

disp(strcat("To run protocol 2 centred on the same position, but with the inverse contrast use the 'peak_frame' value = ", string(inv_index)))

end 