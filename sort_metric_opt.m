function [s1, metric1] = sort_metric_opt(A),

A = transpose(A);

n = size(A,1);

p_vec = perms(flip(1:n));

metric1 = 1e7;

for i=1:size(p_vec,1),

    m = sum(diag(A(p_vec(i,:),1:n)));

    if m < metric1,

        metric1 = m;

        s1 = p_vec(i,:);

    end

end

return