function [ref_image] = l1smoothing(image, splabel, param)
% [ref_image] = l1smoothing(image, splabel)
% Usage: smoothing input image with L1 optimization
% Input:
%   - image: input image
%   - splabel: super-pixel label for each pixel
%   - param:
%       .alpha       [20] local sparseness weight
%       .beta        [0.01] global sparseness weight
%       .theta       [50] image approximation term
%       .lambda      [120] regularization term weight
%       .itr_num     [4] number of iterations in optimization
%       .window_size [5] size of local sparseness windows
% Output:
%   - ref_image: output image after smoothing
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

width = size(image, 2); height = size(image, 1); pixel_num = width * height;
image = double(image);
r = image(:,:,1); r = r(:);
g = image(:,:,2); g = g(:);
b = image(:,:,3); b = b(:);

% construct first term matrix: local sparse
fprintf('Construct local sparse matrix...\n');
A = windowvar(image, window_size, param); A = alpha * A;

% construct second term matrix: global sparse
fprintf('Construct global sparse matrix...\n');
B = spvar(image, splabel);
B = beta * B;

target = [r; g; b];

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
end

ref_image = reshape(ref, height, width, 3);
ref_image = uint8(ref_image);
end


function d = shrink(v, lambda)
index = find(v < 0);
d = max(abs(v) - lambda, 0);
d(index) = -1 * d(index);
end




