function [M] = windowvariation(image, window_size, param)
% [M] = windowvariation(image, window_size)
% Usage: generate a matrix used to get the total varation for each pixel 
%	in a window centered at it.
% Input:
%	- image: original image
%	- window_size: half size of the window
%	- param: parameters
% Ouput:
% 	- M: generated matrix (sparse)
if isfield(param, 'mu')
	mu = param.mu;
else
	mu = 10.0;
end

if isfield(param, 'ga')
	ga = param.ga;
else
	ga = 120.0;
end

if isfield(param, 'sigma')
	sigma = param.sigma;
else
	sigma = 0.5;
end

cform = makecform('srgb2lab');
image_lab = applycform(uint8(image),cform);
image_lab = double(image_lab);

% image = double(image);
% image = image / 255.0;

height = size(image, 1); 
width = size(image,2);
pixel_num = height * width;

est_pair_num = pixel_num * window_size * window_size;
row = zeros(2*est_pair_num, 1);
col = zeros(2*est_pair_num, 1);
val = zeros(2*est_pair_num, 1);

count = 1;
row_count = 1;
small = @(x,y) min([x,y]);

chrom = image_lab(:,:,1) / 100.0 ;
chrom_r = image_lab(:,:,2) / 220.0;
chrom_g = image_lab(:,:,3) / 220.0;

chrom = mu * chrom; chrom = chrom(:); % 10: best 
chrom_r = ga * chrom_r; chrom_r = chrom_r(:);
chrom_g = ga * chrom_g; chrom_g = chrom_g(:);


for i = 1 : width
	for j = 1 : height
		p = i : small(width, i+window_size);
		q = j : small(height, j + window_size);
		if(isempty(p) || isempty(q)) continue; end

		[x,y] = meshgrid(p,q);
		x = x(:); y = y(:);
		if(length(x) < 2) continue; end

		center_id = (i - 1) * height + j;
		neigh_id = (x - 1) * height + y;
		neigh_id = neigh_id(2:end);
		pair_num = length(neigh_id);
		index = count:(count+pair_num-1);
		row_index = row_count : (row_count + pair_num - 1);
		row(index) = row_index;
		col(index) = center_id;
		index_1 = (count+pair_num) : (count + 2 * pair_num -1);
		row(index_1) = row_index;
		col(index_1) = neigh_id;

		color_diff = [chrom(neigh_id) - chrom(center_id) ...
					  chrom_r(neigh_id) - chrom_r(center_id) ...
					  chrom_g(neigh_id) - chrom_g(center_id)];
		color_diff = sum(color_diff.^2, 2);
		val(index) = abs(color_diff);
		val(index_1) = -1.0 * abs(color_diff);

		count = count + 2 * pair_num;
		row_count = row_count + pair_num; 
	end
end

row = row(1:count-1);
col = col(1:count-1);
val = val(1:count-1);

% calculate variance 
index_1 = find(val > 0);
val(index_1) = exp(-0.5 * val(index_1)); 
index_1 = find(val < 0);
val(index_1) = -1.0 * exp(sigma * val(index_1));


row_1 = row + length(row) / 2;
col_1 = col + pixel_num ;

row_2 = row_1 + length(row) / 2;
col_2 = col_1 + pixel_num;

final_row = [row;row_1;row_2];
final_col = [col;col_1;col_2];
final_val = [val;val;val];

M = sparse(final_row, final_col, final_val);

