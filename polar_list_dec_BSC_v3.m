%This function implements the list decoder of a polar code
%inputs are as follows:
%LV_r: beliefs received via the BSC channel. Length of this belief is the
%codeword length
% K: Number of information bits
% Rseq: Reliability sequence
% L: Number of surviving paths in the list
% emd_mat: a matrix showing the sequence of edges and move directions in order 
%of their appearance in the decoding process. This sequence is constructed
%once using the function buld_emd and is used throughout the decoding
%hence build_emd(N) must be run before this decoder function
%N is the codeword length
%eb_loc shown the position of each message in the path(s) and is also constructed
%by running build_emd(N) prior to running this decoder

%function outputs are as follows:
%x_mat: a matrix containing L sequences of K bit decoded messages
%metric_vec: metrics corresponding to each path

%Note that this list decoder does not use CRC by itself. One needs to
%implement CRC separately and add it to this scheme

function [x_mat, metric_vec] = polar_list_dec_BSC_v3(Lv_r,K,Rseq,L,emd_mat,eb_loc),

%Codeword length
N = length(Lv_r);

%find the message length for each edge
de = floor(log2(1:2*N-1));
MLE = N ./ (2.^de);

%mark the edges corresponding to the K information bits
IIB = zeros(1,2*N-1);
IIB(N-1+Rseq(N-K+1:N)) = 1;

%initiate surviving paths and metrics
path_mat = Lv_r;
metric_vec = [0];

total_moves = size(emd_mat,2);

edge_vec = emd_mat(1,:);

md_vec = emd_mat(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for mov = 1:total_moves,

    CE = edge_vec(mov); CMD = md_vec(mov);

    if (CE >= N && CMD == 2 && IIB(CE) == 0) %moving upward from a frozen bit
        
        %fetch the downward beliefs (DWBs) over this edge through all paths. These beliefs are
        %located at the column of the current path_mat 

        DWB_vec = path_mat(:,end);

        %Generate the LLRs for these beliefs (they are currently L-values between 0 and 1)
        %to ensure that we do not get INF and NAN in the logarithm, we
        %first clip DWB:
        eps = 1e-4;
        DWB_vec = min(eps,max(1-eps,DWB_vec));
        L_vec = log(DWB_vec ./ (1-DWB_vec));

        %we extend all paths by adding an upward belief which is 0 for the
        %frozen bit
        NP = size(path_mat,1);
        path_mat = [path_mat, 0*ones(NP,1)];

        %We penalize a path only if the DWB is biased toward 1 
        % (i.e. DWB > 0.5)
        MP_vec = abs(L_vec) .* (DWB_vec > 0.5); %max(L_vec,0);
        metric_vec = metric_vec + MP_vec;

        %End of this case

    elseif (CE >= N && CMD == 2 && IIB(CE) == 1) %moving upward from an information bit

        %fetch the downward beliefs (DWBs) over this edge through all paths. These beliefs are
        %located at the column of the current path_mat 

        DWB_vec = path_mat(:,end);

        %Generate the LLRs for these beliefs (they are currently L-values between 0 and 1)
        %to ensure that we do not get INF and NAN in the logarithm, we
        %first clip DWB:
        eps = 1e-4;
        DWB_vec = max(eps,min(1-eps,DWB_vec));
        L_vec = log(DWB_vec ./ (1-DWB_vec));

        %we extend all paths by adding both possible upward beliefs (0 and 1) 
        NP = size(path_mat,1);
        path_ext = [0*ones(NP,1) ; 1*ones(NP,1)];
        
        path_mat = [path_mat; path_mat];
        path_mat = [path_mat, path_ext];

        %We penalize a path only if the DWB is biased against its extension 
        % (i.e. DWB > 0.5 for extension 0 and DWB < 0.5 for extension 1)
        MP_vec1 = abs(L_vec) .* (DWB_vec > 0.5);
        MP_vec2 = abs(L_vec) .* (DWB_vec < 0.5);
        MP_vec = [MP_vec1; MP_vec2];

        metric_vec = [metric_vec; metric_vec];
        metric_vec = metric_vec + MP_vec;

        %End of this case


    elseif (CE < N && CMD == 2) %moving upward from a non-leaf edge
        %the serving edges are left child edge (LCE) and right child edge RCE
        %first we need to extract the upward beliefs of LCE and RCE from
        %the paths
        LCE = 2*CE;
        RCE = 2*CE+1;

        loc1 = eb_loc(2,LCE);
        loc2 = eb_loc(2,RCE);
        
        %message length of LCE and RCE
        mln = MLE(LCE);
        
        
        %fetch the upward messages of the LCE from all paths
        UMLCE = path_mat(:,loc1:loc1+mln-1);

        %fetch the upward messages of the RCE from all paths
        UMRCE = path_mat(:,loc2:loc2+mln-1);

        %Now combine the messages to find UWB for CE
        path_ext = [mod(UMLCE+UMRCE,2), UMRCE];
        
        path_mat = [path_mat, path_ext];

        %End of this case


    elseif (mod(CE,2) == 0 && CMD == 1) %moving downward on a left child edge
        %The serving edge is the parent edge
        PE = floor(CE/2);

        %fetch the DWB of the PE
        loc = eb_loc(1,PE);
        mln = MLE(PE);
        DWMPE = path_mat(:,loc:loc+mln-1);

        %Split the downward belief of the PE
        DLP = DWMPE(:,1:mln/2);
        DRP = DWMPE(:,mln/2+1:end);
        
        %Update the belief and pass
        f = DLP .* (1-DRP) + DRP .* (1-DLP);

        path_mat = [path_mat, f];

        %end of this case


    elseif (mod(CE,2) == 1 && CMD == 1) %moving downward on a right child edge

        %the serving edges are the PE and the left sibling edge (LSE)
        PE = floor(CE/2);
        LSE = CE-1;

        %fetch the DWB of the PE
        loc = eb_loc(1,PE);
        mln = MLE(PE);
        DWMPE = path_mat(:,loc:loc+mln-1);

        %Split the downward belief of the PE
        DLP = DWMPE(:,1:mln/2);
        DRP = DWMPE(:,mln/2+1:end);

        %fetch the UWB of the LSE
        loc = eb_loc(2,LSE);
        mln = MLE(LSE);

        ULS = path_mat(:,loc:loc+mln-1);
        
        %calculate beliefs and pass
        fup = ULS .* (1-DLP) + (1-ULS) .* DLP;

        g = (fup .* DRP)./(fup .* DRP + (1-fup) .* (1-DRP)) ;

        path_mat = [path_mat, g];

        %end of this case

       
    else
        disp('unexcpected case in list decoding process');
    end

    %If the number of surviving paths exceeds L, reduce them to L
    % Note that paths with smaller metrics are preferred, since metrics
    % are increased by adding penalties

    if size(path_mat,1) > L,

        
        [i, j] = sort(metric_vec);

        path_mat = path_mat(j(1:L),:);

        metric_vec = metric_vec(j(1:L));

        
    end

end

%Pick the message bits corresponding to each path

ji = Rseq(N-K+1:N)+N-1;

m_loc = eb_loc(2,ji);

x_mat = path_mat(:,m_loc);


return

