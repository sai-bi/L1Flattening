image_path = 'C:\Users\bisai\Documents\research\dataset\l1flttening\53.png';
image = imread(image_path);
sp_path = 'C:\Users\bisai\Documents\Visual Studio 2013\Projects\L1Flattening\x64\Release\L1Flattening.exe';
sp_output_path = 'C:\Users\bisai\Documents\research\L1\matlab\sp.txt';
command = [sp_path ' -i ' image_path ' -o ' sp_output_path];
system(command);
splable = csvread(sp_output_path);