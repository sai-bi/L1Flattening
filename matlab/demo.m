addpath('./../src');
image_path = '..\data\2.png';
image = imread(image_path);
% use default parameters for image flattening
param = struct(); 
flat_image = l1flattening(image, param);
imwrite(flat_image, '..\data\2-flat.png');

% use default parameters for edge-preserving smoothing
param.local_param.edge_preserving = true; 
flat_image = l1flattening(image, param);
imwrite(flat_image, '..\data\2-smooth.png');
