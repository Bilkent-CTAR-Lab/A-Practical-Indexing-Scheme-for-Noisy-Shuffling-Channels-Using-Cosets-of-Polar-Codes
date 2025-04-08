clc;
clear;

logtanh = log(tanh(0.001:0.001:4.001));

atanhexp = atanh(exp(0:-0.001:-8.001));



%Coset types: (1:explicit) (2:random) (3:designed) 
% (4:random in search area)

coset_type = 2;

pd_vec = [];


%Parameters of the outer code: RS code
m = 8; %8;           % Number of bits per symbol
n = 2^m - 1;     % Codeword length
%k = 225; %225;           % Message length

List_size = 1;

%Number of information bits per whole block
%lbits = m*k;

%number of segments
seg_num = 16;

%number of bits required to encode the index of each segment
ind_bits_num = ceil(log2(seg_num));

%number of information bits in each slice
% k_seg = ceil(n*m/seg_num);
k_seg = ceil(n*m/seg_num);

%Outer code : polar codeword length
n_seg = 256; %128;

dv = 4; dc = 8; 

%[H_orig, G_sys] = build_incomplete_Gallager(dv,dc,n_seg,k_seg);


%[emd_mat,eb_loc] = build_emd(n_seg);

%Load the 5G standard reliability sequence for 1024
%rel_seq1024;


%Extract reliability sequence for n_seg
%[ir,jr] = find(Rseq1024 <= n_seg);

%Rseq = Rseq1024(jr);

%Build polar transform matrix
%n_polar = floor(log2(n_seg));

%G = build_G(n_polar);

%This flag bit denotes whether explicit indexing is used
%ind = 0

k_vec = 225:10:225;

p_bsc_vec = 0.035:0.005:0.050;

BER0 = ones(length(k_vec),length(p_bsc_vec));

BER1 = ones(length(k_vec),length(p_bsc_vec));

FER0 = ones(length(k_vec),length(p_bsc_vec));

FER1 = ones(length(k_vec),length(p_bsc_vec));

SER0 = ones(length(k_vec),length(p_bsc_vec));

SER1 = ones(length(k_vec),length(p_bsc_vec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metric_mat = 1e5*ones(seg_num,seg_num);

if(coset_type == 3),

    %load good_cosets_m6;
    
    %cn_e_vec = co_vec_6;

    %nfb = size(cn_e_vec,2);

    %coset_mat = e_vec_6;

    
    %load e_vec_63_128.mat;
    
    %coset_mat = e_vec_6;

    %load good_example_32_4cosets;

    %coset_mat = e_vec;

    
    %load e_vec_128_64_32.mat;

    %coset_mat = e_vec_32;

    %load e_vec_good_42_64_9_12fb.mat;

    %coset_mat = e_vec_42_64_9;

    %load good_cosets_m6_nf6;

    %coset_mat = e6;

    load e_vec_Hamming_64_42_9_nf12;

    coset_mat = e_Hamming;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sent_blocks = 2000;

%load test_data;

%load u_mat_load;

%load z_load;

symb_errors_mat = zeros(length(p_bsc_vec),sent_blocks);

for ind = 0:0,

    %AWGN variance
    %SNR_db = 1;

    for ik = 1:length(k_vec),

        k = k_vec(ik);

        lbits = m*k;

        % udl = floor(length(u_data)/lbits);
        % 
        % u_data_vec = reshape(u_data(1:lbits*udl),udl,lbits);
        % 
        % zdl = floor(length(z_data)/n_seg);
        % 
        % z_data_vec = reshape(z_data(1:zdl*n_seg),zdl,n_seg);

        

        %u_mat = u_data_vec(1:sent_blocks,:); 
        
        u_mat = round(rand(sent_blocks,lbits));

        for ipb = 1:length(p_bsc_vec),

            p_bsc = p_bsc_vec(ipb);

            %sigma_noise = sqrt(0.5*10^(-SNR_db/10));

            bit_errors = 0;

            frame_errors = 0;

            symbol_errors = 0;

            nde = 0;

            

            %coset_mat = rand(seg_num , n_seg) < 0.1;
            
            %Hcs = round(hadamard(n_seg) > 0);

            %coset_mat = Hcs(2:seg_num+1,:);

            %ercs = ones(seg_num,1) * round(rand(1,n_seg));

            %coset_mat = rem(Hcs(1:seg_num,:) + ercs, 2);

            %coset_mat = seqs_max_dist(seg_num,n_seg,10,20,0.1);

            %coset_mat = select_coset_leaders(seg_num,80,Rseq);

            
            for sb = 1:sent_blocks,

                [H_orig, G_sys] = build_incomplete_Gallager(dv,dc,n_seg,k_seg);

                %block = sb

                %Define interleavers and de-interleavers for matched decoder scheme
                %[intliv_mat,deintliv_mat] = define_interleavers(seg_num,n_seg);

                %Define coset leaders for matched decoder scheme

                %%%%%%%%%%%%%%%% Explicit indexing
                
                if(coset_type == 1),

                    coset_mat = [];
    
                    for snco = 0:seg_num-1,
                        
                        u_snco = zeros(1,n_seg);
    
                        selected_pos = Rseq(n_seg-k_seg-ind_bits_num+1:n_seg-k_seg);
    
                        u_snco(selected_pos) = gf2bit(snco,ind_bits_num);
    
                        selected_pos = Rseq(1:n_seg-k_seg); %-ind_bits_num+1:n_seg-k_seg);
    
                        % 
                        cs_leader = rem(u_snco * G , 2);
    
                        coset_mat = [coset_mat ; cs_leader];
    
                    end

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Random Cosets

                if(coset_type == 2),
                    coset_mat = round(rand(seg_num,n_seg) < 0.5);
                    %coset_mat = zeros(seg_num,n_seg);
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%% Random in search area
                %load dms_app_63_128_16.mat e_vec;

                if(coset_type == 4),

                    %naaa = size(e_vec,1);

                    %coset_mat = e_vec(randperm(naaa, seg_num),:);

                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

                %generating information bits
                %u = u_mat(sb,:); 
                u = round(rand(1,lbits) < 0.5);

                %Converting bits to m-bit symbols
                unb = bit2gf(u,m);

                %converting symbols to MATLAB's GF(2^m) data format
                msg = gf(unb,m);

                %Encoding by RS code
                code = rsenc(msg,n,k);

                %Converting the codeword from MATLAB's GF(2^m) format to integer format
                codeint = zeros(1,n);

                for j=0:2^m-1,

                    codeint = codeint + j*(code == j);

                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Converting the symbols sequences to bit sequence
                x = gf2bit(codeint,m);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Slicing the sequence into segments

                %padding the input sequence to get an integer number of segments
                if length(x) < k_seg*seg_num,

                    x = [x, zeros(1,k_seg*seg_num-length(x))];
                end

                %encoding each slice by a polar code

                ldpc_cws = [];

                for si = 1:seg_num,

                    xi = x(1+(si-1)*k_seg:si*k_seg);

                    if ind == 1,

                        xi = [xi, gf2bit(si-1,ind_bits_num)];
                    end

                    %Calling the polar code function
                    xi_out = rem(xi*G_sys, 2);%polar_encode(xi,n_seg,G,Rseq);
                    %%%%%

                    if ind == 1,
                        %polar_cws = [polar_cws; xi_out];
                    end


                    %interleave each segment if matched decoding is employed
                    if ind == 0,

                        %polar_cws = [polar_cws; xi_out(intliv_mat(si,:))];

                        ldpc_cws = [ldpc_cws; rem(xi_out+coset_mat(si,:),2)];
                    end


                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Passing through a BSC channel
                
                %z_mat = round(z_data_vec((sb-1)*seg_num+1:sb*seg_num,:) < p_bsc); 
                
                z_mat = (rand(size(ldpc_cws)) < p_bsc);

                y_mat = rem(ldpc_cws+z_mat, 2);

                %Shuffle the noisy observations
                shuf = randperm(seg_num);

                y_mat = y_mat(shuf,:);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Decoding the segments and detecting the indexes to put them in order
                
                x_rec_mat = zeros(seg_num,k_seg);

                %decoding in case of explicit indexing
                if ind == 1,
                    % for si=1:seg_num,
                    % 
                    %     %calling the polar decoder
                    %     k_val = k_seg+ind_bits_num;
                    % 
                    %     %a priori probability of message bits being "1"
                    %     %f_vec = 0.5*ones(1,n_seg);
                    %     %a priori reliability of frozen bits being "1" is zero
                    %     %f_vec(Rseq(1:n_seg-k_seg)) = 0;
                    % 
                    %     out_rel_vec = (1-p_bsc)*y_mat(si,:)+(p_bsc)*(1-y_mat(si,:));
                    % 
                    %     %xi_rec = polar_dec_BSC(y_mat(si,:),k_seg+ind_bits_num,Rseq);
                    % 
                    %     xi_rec = polar_dec_BSC(out_rel_vec,k_val,Rseq);
                    % 
                    %     %[out_rec, xi_rec] = polar_decode(out_rel_vec,f_vec);
                    % 
                    %     xi_ind_bits = xi_rec(end-ind_bits_num+1:end);
                    % 
                    %     xi_index = bit2gf(xi_ind_bits,ind_bits_num) + 1;
                    % 
                    %     x_rec_mat(xi_index,:) = xi_rec(1:end-ind_bits_num);
                    % 
                    %     %x_rec_mat(xi_index,:) = xi_rec(Rseq(n_seg-k_seg+1:end));
                    % 
                    % end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %decoding in case of matched decoder detection
                if ind == 0;
                    %keep track of indexes that are already detected
                    detected_inds = zeros(1,seg_num);

                    vd1_all = [];

                    rr1_all = [];

                    for si = 1:seg_num,

                        
                        for sj = 1:seg_num,

                            if 1,

                                %calling the polar decoder
                                
                                y_deliv = rem(y_mat(si,:)+coset_mat(sj,:),2);

                                %out_rel_vec = (1-p_bsc)*y_deliv + (p_bsc)*(1-y_deliv);

                                out_rel_vec = log((1-p_bsc)/p_bsc)*(1-2*y_deliv);

                                %[vd1_mat, metric1_vec] = polar_list_dec_BSC_v3(out_rel_vec,k_seg,Rseq,List_size,emd_mat,eb_loc);
                                
                                ldpc_iters = 10;

                                [xd1_mat,synd] = ldpc_dec_sp(out_rel_vec, H_orig, ldpc_iters);

                                                                
                                %[xd1_mat,synd] = ldpc_dec_sp_lookup(out_rel_vec, H_orig, ldpc_iters,logtanh,atanhexp);
                                                              
                                %xd1_mat = zeros(List_size,n_seg);
                
                                %xd1_mat(:,Rseq(n_seg-k_seg+1:n_seg)) = vd1_mat;
                
                                vv1_mat = rem(xd1_mat*G_sys,2);
                
                                dd1_vec = sum(transpose(rem(ones(List_size,1)*y_deliv + vv1_mat , 2)));
                
                                [dd1, loc_dd1] = min(dd1_vec);
                
                                vd1_all = [vd1_all; xd1_mat(loc_dd1,:)];
                                
                                %input('press a key');

                                vv1_best = vv1_mat(loc_dd1,:);

                                rr1_best = rem(vv1_best + coset_mat(sj,:),2);

                                rr1_all = [rr1_all; rr1_best]; 

                                metric_mat(si, sj) = dd1;


                            end

                        end

                        
                    end

                    [sort_out, metric_out] = sort_metric(metric_mat);

                    %sort_out = shuf

                    %ss = [sort_out;shuf]

                    %input('press a key');

                    %[sort_out, metric_out] = sort_metric_opt(metric_mat);

                    if(~isequal(sort_out,shuf)),

                        nde = nde + 1

                    end

                    %ss = [sort_out;shuf]

                    %input('press a key');

                    
                    %x_rec_mat = [];

                    rr1_dec = [];

                    for iii = 1:seg_num,

                        xxx = vd1_all((iii-1)*seg_num + sort_out(iii), :);

                        %x_rec_mat = [x_rec_mat; xxx];

                        x_rec_mat(sort_out(iii),:) = xxx;

                        rrr = rr1_all((iii-1)*seg_num + sort_out(iii), :);

                        rr1_dec = [rr1_dec; rrr];

                    end


                end

                

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %reshape polar decoders ordered output
                x_vec = reshape(transpose(x_rec_mat),1,prod(size(x_rec_mat)));

                %ldpc_ber = sum(x ~= x_vec)/length(x_vec)

                %input('press a key');

                %Remove the padding bits that were added at the encoder for segmentation
                x_vec = x_vec(1:n*m);

%                 if(isequal(shuf,sort_out)),
% 
%                     input('press a key');
% 
%                 end

                %change the received bits from the polar decoder into symbols
                x_vec_symbols = bit2gf(x_vec,m);

                %Remove the padding bits that were added at the encoder for segmentation
                %x_vec_unpad = x_vec_symbols(1:n);

                %change the format of symbols into MATLAB's GF(2^m) symbols
                code_rec = gf(x_vec_symbols,m);

                %decode using RS decoder
                [rxcode,cnumerr] = rsdec(code_rec,n,k);

                symb_errors_mat(ipb,sb) = sum(code_rec ~= code);

                %input('press a key');

                rxint = zeros(1,k);

                for j=0:2^m-1,

                    rxint = rxint + j*(rxcode == j);

                end


                udec = gf2bit(rxint,m);

                %[u;udec]

                bit_errors = bit_errors + sum(u ~= udec)

                bx = sum(u ~= udec);

                if(bx > 0),
                    frame_errors = frame_errors+1;
                end

                symbol_errors = symbol_errors + sum(unb ~= rxint);

                %sb_ber = [sb, bit_errors/sb/lbits]

                %sb_fer = [sb, frame_errors]

                sb_ser = [sb, symbol_errors];

                if(rem(sb,10) == 0),
                    
                    if(rem(sb,20)==0),
                        toc;
                    else
                        tic
                    end

                    sb_ind_k_p_fe = [sb, ind, k, p_bsc, frame_errors]
                    %input('press a key');
                end

            end

            BER = bit_errors / lbits / sent_blocks;

            FER = frame_errors / sent_blocks;

            SER = symbol_errors / sent_blocks / length(rxint);

            if ind == 0,

                BER0(ik,ipb) = BER;
                FER0(ik,ipb) = FER;
                SER0(ik,ipb) = SER;

            end

            if ind == 1,

                BER1(ik,ipb) = BER;
                FER1(ik,ipb) = FER;
                SER1(ik,ipb) = SER;
                
            end

            pd_vec = [pd_vec, nde/sent_blocks]


        end

    end

end


save BSC_ldpc_results1.mat;


