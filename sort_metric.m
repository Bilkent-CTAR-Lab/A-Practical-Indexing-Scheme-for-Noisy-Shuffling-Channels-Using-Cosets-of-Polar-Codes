function [s1, metric1] = sort_metric(A),

metric1 = 0;

A = transpose(A);

m = max(max(A)) + 1;

s = zeros(1,size(A,2));

for i = 1:size(A,2),

    [B, J] = sort(A);

    %B = B

    d = B(2,:) - B(1,:);

    [val, loc] = sort(d);

    for kk = size(A,2):-1:1,

        jj = J(1,loc(kk));

        loc1 = loc(kk);

        if(s(loc1) == 0),

            s(loc1) = jj;

            metric1 = metric1 + B(1,loc1);

            break;
        end

    end

    %s(loc) = jj;

    %metric1 = metric1 + B(1,loc);

    A(jj,:) = m;

    A(:,loc1) = m;

    %input('press a key');

end

s1 = s;

return