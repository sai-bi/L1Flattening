function [ref_image] = l1smoothing(image, splabel, param)
% [ref_image] = l1smoothing(image, splabel)
% Usage: smoothing input image with L1 optimization
% Input:
%   - image: input image
%   - splabel: super-pixel label for each pixel
%   - param: parameters, alpha, beta, theta, lambda
% Output:
%   - ref_image: output image after smoothing

width = size(image, 2);
height = size(image, 1);
pixel_num = width * height;
image = double(image);
r = image(:,:,1); r = r(:);
g = image(:,:,2); g = g(:);
b = image(:,:,3); b = b(:);

% construct first term matrix: local sparse
fprintf('Construct local sparse matrix...\n');
if ~exist('ori_A', 'var')
    window_size = 5;
    ori_A = windowvariation(image, window_size, param);
end

% construct second term matrix: global sparse
fprintf('Construct global sparse matrix...\n');
B = spvariation(image, splabel);
target = [r; g; b];

% alpha = 20; beta = 0.01; theta = 50; lambda = 120;
% check parameters
if isfield(param, 'alpha')
    alpha = param.alpha;
else
    alpha = 20;
end

if isfield(param, 'beta')
    beta = param.beta;
else
    beta = 0.01;
end

if isfield(param, 'theta')
    theta = param.theta;
else
    theta = 50;
end

if isfield(param, 'lambda')
    lambda = param.lambda;
else
    lambda = 120;
end

A = alpha * ori_A; clear ori_A;
B = beta * B;

fprintf('Calculate left hand matrix...\n');
left_hand = lambda * (A' * A) + (B' * B) + ...
        theta * sparse(1:3*pixel_num, 1:3*pixel_num, ones(1,3*pixel_num));

itr_num = 4; 

ref = zeros(pixel_num*3,1);
old_ref = target;
threshold = 0.001;
b_1 = zeros(size(A,1),1);
d_1 = zeros(size(A,1),1);
b_2 = zeros(size(B,1),1);
d_2 = zeros(size(B,1),1);
for i = 1 : itr_num
    fprintf('Iteration %d out of %d...\n',i, itr_num);
    if(norm(ref - old_ref) < threshold)
        break;
    end 
    old_ref = ref;

    right_hand = theta * target + lambda * (A' * (d_1 - b_1) + B' * (d_2 - b_2));

    ref = left_hand \ right_hand;
    temp_1 = A * ref;
    temp_2 = B * ref;
    d_1 = shrink(temp_1 + b_1, 1.0 / lambda);
    d_2 = shrink(temp_2 + b_2, 1.0 / lambda);
    b_1 = b_1 + temp_1 - d_1;
    b_2 = b_2 + temp_2 - d_2;

    % part_1 = sum(abs(temp_1)) / (alpha);
    % part_2 = sum(abs(temp_2)) / (beta);
    % part_3 = sum((ref-target).^2);
    % total = alpha * part_1 + beta * part_2 + theta * part_3;

    % fprintf('Iter %d, total: %f\n',i, total);
    % fprintf('1: %f, 2: %f, 3: %f \n', part_1, part_2, part_3);
end

ref_image = reshape(ref, height, width, 3);