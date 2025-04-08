function [H_orig, G_sys] = build_incomplete_Gallager(dv,dc,n,ki),

l = n/dc;

inv_flags = [];

force_flip = 0;

A = zeros(l, dc*l);

for i = 1:l,
    A(i,dc*i-dc+1:dc*i) = ones(1,dc);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_orig = A;

for j=2:dv,

    H_orig = [H_orig; A(:,randperm(dc*l))];

end

H_orig = H_orig(1:n-ki,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m, n] = size(H_orig);

H_sys = H_orig;

inv_flag = 1;

rni = 0;

for i = 1:m,

    inv_flag_i = 1;

    if(H_sys(i,i) == 0),

        inv_flag_i = 0;

        for t = i+1:m,

            if(H_sys(t,i) == 1),

                inv_flag_i = 1;

                temp_vec = H_sys(i,:);
                H_sys(i,:) = H_sys(t,:);
                H_sys(t,:) = temp_vec;

                temp_vec = H_orig(i,:);
                H_orig(i,:) = H_orig(t,:);
                H_orig(t,:) = temp_vec;
            end

        end

        
    end

    if inv_flag_i == 0,

        force_flip = force_flip + 1;
        
        H_sys(i,i) = 1;

       
        H_orig(i,i) = mod(H_orig(i,i)+1,2);

        
        inv_flag_i = 1;

    end

    if inv_flag_i ~= 0,
    
        for j = 1:size(H_sys,1),

            if j ~= i,
    
                if H_sys(j,i) == 1
                
                    V = H_sys(j,:) + H_sys(i,:);
    
                    H_sys(j,:) = V - 2*floor(V/2);
    
                end
    
            end

        end
    end

    inv_flags = [inv_flags, inv_flag_i];

    rni = rni + inv_flag_i;

    ir = [i, rni];

    inv_flag = inv_flag * inv_flag_i;

   
end


mf = size(H_sys,1);

G_sys = [transpose(H_sys(:,mf+1:n)), eye(n-mf)];


return