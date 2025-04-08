function G = build_G(n),

if n==1,

    G = [1 0;1 1];

else

    G = [build_G(n-1) zeros(2^(n-1));build_G(n-1) build_G(n-1)];

end

return;