image_path = 'C:\Users\bisai\Documents\research\dataset\l1flttening\53.png';
image = imread(image_path);
sp_path = '..\bin\L1Flattening.exe';
sp_output_path = '..\bin\sp.txt';
command = [sp_path ' -i ' image_path ' -o ' sp_output_path];
system(command);
splabel = csvread(sp_output_path);

[ref_image] = l1smoothing(image, splabel, 0);
imwrite(uint8(ref_image), '..\bin\53-flat.png');
