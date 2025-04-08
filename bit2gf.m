function unb = bit2gf(u,m),

L = length(u)/m;

if(ceil(L) ~= L),

    disp('Invalid inputs to bit2gf(u,m): length of u must be a multiple of m');

end

x = reshape(u,m,L);

unb = zeros(1,L);

for i = 0:m-1,

    unb = unb + x(i+1,:)*(2^(m-1-i));

end

return;