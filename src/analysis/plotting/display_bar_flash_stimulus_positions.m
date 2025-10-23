
% Plot the bar flash stimulus:

% Load the bar flash pattern "0011" in the P2 experiments from October 2025
% onwards.
% pat = pattern.Pats(:, :, 2:end);

orient_id = 7;
close 

A = cell(1, 11);

for k = 1:11
    frame_id = k+(11*orient_id);
    f = pat(:, :, frame_id);
    f2 = f;
    f2(f2>6) = 1;
    f2(f2>1) = 0;
    A{k} = f2;
end 

% Build a label image: 0 = background, k = bar from frame k.
L = zeros(size(A{1}), 'uint8');
for k = 1:numel(A)
    L(A{k} ~= 0) = k;   % if bars overlap, the later one (higher k) wins
end

cmap = create_red_blue_cmap();

% Turn labels into an RGB image (background = black).
% (ind2rgb expects indices starting at 1, so shift by +1 and prepend bg color)
RGB = ind2rgb(double(L)+1, [1 1 1; cmap]);

figure
imshow(RGB)
f = gcf;
f.Position = [2088   -165  2392  1094]; %[13  35 1771 937];



