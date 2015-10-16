image_path = '..\data\2.png';
image = imread(image_path);
sp_path = '..\bin\L1Flattening.exe';
sp_output_path = '..\data\sp.txt';
command = [sp_path ' -i ' image_path ' -o ' sp_output_path];
system(command);
splabel = csvread(sp_output_path);

% use default parameters for image flattening
param = struct(); 
flat_image = l1flattening(image, splabel, param);
imwrite(flat_image, '..\data\2-flat.png');

% use default parameters for edge-preserving smoothing
param.local_param.edge_preserving = true; 
flat_image = l1flattening(image, splabel, param);
imwrite(flat_image, '..\data\2-smooth.png');
