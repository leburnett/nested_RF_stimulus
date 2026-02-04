function find_inv_peak_frame(peak_frame)
% FIND_INV_PEAK_FRAME  Calculate peak frame number for opposite contrast.
%
%   FIND_INV_PEAK_FRAME(PEAK_FRAME) computes and displays the frame number
%   needed to run Protocol 2 at the same spatial location but with the
%   opposite contrast (ON->OFF or OFF->ON).
%
%   INPUT:
%     peak_frame - Frame number from Protocol 1 that showed best response
%                  (typically 1-256 for the 6px flash grid)
%
%   FRAME NUMBERING:
%     In Protocol 1, frames 1-128 correspond to one contrast (e.g., ON),
%     and frames 129-256 correspond to the opposite contrast (OFF) at
%     the same spatial locations. This function maps between them.
%
%   USAGE:
%     After running Protocol 2 with one contrast preference, use this
%     function to find the frame number for testing the same neuron's
%     response to the opposite contrast at the identical location.
%
%   EXAMPLE:
%     % If ON flash at frame 45 gave best response:
%     find_inv_peak_frame(45)
%     % Displays: "To run protocol 2... use the 'peak_frame' value = 173"
%
%   See also GENERATE_PROTOCOL2, PATT_FRAME_TO_COORD

if peak_frame < 130
    inv_index = peak_frame+128;
else
    inv_index = peak_frame-128;
end 

disp(strcat("To run protocol 2 centred on the same position, but with the inverse contrast use the 'peak_frame' value = ", string(inv_index)))

end 