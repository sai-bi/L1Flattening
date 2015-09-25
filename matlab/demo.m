image_path = '..\data\1.png';
image = imread(image_path);
sp_path = '..\bin\L1Flattening.exe';
sp_output_path = '..\data\sp.txt';
command = [sp_path ' -i ' image_path ' -o ' sp_output_path];
system(command);
splabel = csvread(sp_output_path);

param = struct(); % use default parameters
flat_image = l1flattening(image, splabel, param);
imwrite(flat_image, '..\data\1-flat.png');
figure;
imshow(flat_image);
