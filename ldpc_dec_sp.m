function [u_dec,synd] = ldpc_dec_sp(Li, H_orig, Max_iter),

Large_LLR = 1e10;

[m, n] = size(H_orig);

mv = H_orig .* (ones(m,1)*Li);

mvc = mv;



for t = 1:Max_iter,

    %Check node to variable node message passing

    mcv = zeros(size(H_orig));

    for cn = 1:m,

        [dum_vec, j_vec] = find(H_orig(cn,:));

        for l = 1:length(j_vec),

            loc = j_vec(l);

            prod1 = 1;

            for lt = 1:length(j_vec),

                if lt ~= l,

                    loct = j_vec(lt);

                    prod1 = prod1 * tanh(mvc(cn,loct)/2);

                 end

            end

            if abs(prod1) ~= 1,
                
                mcv(cn,loc) = log((1+prod1)/(1-prod1));
            else

                mcv(cn,loc) = prod1 * Large_LLR;
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Variable node to check node message passing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mvc = zeros(m,n);

    for vn = 1:n,

        sv = sum(mcv(:,vn));

        [i_vec, dum_vec] = find(H_orig(:,vn));

        for l=1:length(i_vec),

            loc = i_vec(l);

            mvc(loc,vn) = sv - mcv(loc,vn);

        end

    end

    mvc = mvc + mv;
    %%%%

    temp_LLRs = sum(mcv) + Li;

    x_hat_temp = round(temp_LLRs < 0);

    synd_temp = mod(H_orig*transpose(x_hat_temp),2);

    if(sum(synd_temp) == 0),
        %check = 1
        t = Max_iter + 1;
    end


    %t = t

    %max_mcv = max(max(abs(mcv)))

    %max_mvc = max(max(abs(mvc)))

    %input('press a key');
    
end

%input('press a key');

final_LLRs = sum(mcv) + Li;

x_hat = round(final_LLRs < 0);

synd = mod(H_orig*transpose(x_hat),2);

u_dec = x_hat(m+1:n);

return

