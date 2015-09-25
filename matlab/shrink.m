function d = shrink(v, lambda)

% d = zeros(size(v));
% for i = 1 : size(v, 1)
%     if v(i) > lambda
%         d(i) = v(i) - lambda;
%     elseif v(i) < -lambda
%         d(i) = v(i) + lambda;
%     end
% end
index = find(v < 0);
d = max(abs(v) - lambda, 0);
d(index) = -1 * d(index);