function u = gf2bit(ub,m),

ub_temp = ub;

u = [];

for i=1:m,

    u = [rem(ub_temp,2);u];

    ub_temp = floor(ub_temp/2);

end

u = reshape(u,prod(size(u)),1);

u = transpose(u);

return