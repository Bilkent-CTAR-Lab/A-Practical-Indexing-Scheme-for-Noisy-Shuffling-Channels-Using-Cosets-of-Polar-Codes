%this function aims to build the path of message flow in the decoding tree 
% of a polar code. 
% 
% Let us define the following abbreviations:
% DW : Downward
% UW : Upward
% DWM: Downward message
% UWM: Upward message
% CE: Current edge (the edge we are visiting right now)
% NE: Next edge (the edge to be visited next)
% CMD: Current Message Direction (Direction of the current message) 
% NMD: Next Message Direction (Direction of the next message)
% PE : Parent Edge (The edge connected to the CE from the top)
% LCE: Left Child Edge (The edge connected to bottom left of CE)
% RCE: Right Child Edge (The edge connected to bottom right of CE)
%LIFE: Leaf Edge (An edge connected to a leaf node of the tree)
%NLIF: Non-leaf edge

% The tree consists of N+1 edges and N+1 nodes, where N is
% the codeword length. The tree has n=log2(N)+1 levels where level 1 consists
% of 1 node (root node) and level i (i between 2 and n), consists of 2*(i-1)
% nodes. The base of tree (level n) consists of N nodes representing the N
% bits of the codeword. 

%Edges and nodes are labeled as follows. The root node is labeled as node 1 
%left child of node i (i between 1 and N-1) is labeled as node 2*i and right
%child of node i is labeled as node 2*i+1. 

%The edge connected to node numbered i from the top (i between 1 and 2*N-1)
%is labeled as edge numbered i. So note that for all i between 2 and 2*N-1
%edge numbered i connects parent node numbered floor(i/2) to its child node
%numbered i. For i=1, edge numbered 1 is just connected to the top of the
%root node and is responsible for delivering the beliefs (L-values) from
%the channel to the decoding tree (in downward direction) and deliver the
%codeword decoded by the tree from node 1 (in upward direction)

%The aim of running this function is to determine the correct sequence of
%edges visited during the decoding process and the direction of messages
%passed during the decoding process; in attempt to avoid repeating
%calculation of this sequence for every decoding.

%the function also calculates the location of each message in the defined 
%decoding path. Each path in the decoder keeps a complete record of all
%message that are passed through the edges up to this point. In a list
%decoder, when a path is split, its current history is copied to its child
%paths and each of the 2 child paths is extended by its own upward belief,
%corresponding to the decision made regarding the examined information bit.

%The location of each message is stored in the array eb_loc, where
%eb_loc(1,i) denotes the start location in the path where the downward message 
% passed through edge i is stored, and eb_loc(2,i) denotes the start location in the 
% path where the upward message passed through edge i is stored.
%length of messages passed through edge i is l(i) = N/(2^s(i)) where s(i) is the level number 
%in the tree where node i is located and is found as s(i) = floor(log2(i))
%therefore, the message is located from the starting point (denoted by the
%array eb_loc) through l(i)-1 next locations.

%the emd_vec array determines the sequence of edges visited during the
%decoding process and direction of messages passed through them
%the initial point is always a downward message through edge 1 that
%deliveres channel beliefs to the tree. We do not store this move in the
%emd_vec array, since it is facilitates generalizing the decoding algorithm later

%e denotes the edge and md denotes the message direction. There are 4 different 
%cases to consider:

%case 1: CE is a LIFE (CE >= N) and CMD = 1
%Since we have reached a leaf node, next we go up through the
% same edge, i.e. NE = CE and NMD = 2
%
%case 2: CE is not a LIFE (CE < N) and CMD = 1
%then continue downward through the LCE
%i.e. NE = 2*CE, NMD = 1

%case 3: CMD = 2 and CE is the LCE of a PE (i.e. CE = 2*i for i=1:N-1)
%then next move is to pass a DWM through the RCE of the PE, i.e.:
% NE = CE + 1, NMD = 1

%case 4: CMD = 2 and CE is the RCE of a PE (i.e. CE = 2*i+1 for i=1:N-1)
%then next move is to pass a UWM through the PE, i.e.:
% NE = floor(CE/2), NMD = 2

%Case 5: CMD = 2, CE = 1, then the codeword is already delivered through
%the upward message passed by edge numbere 1. Terminate the process

%at each move, we also update the eb_loc matrix by writing down the
%start point in the path for the message being passed in the tree

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%N is the codeword length that is an integer power of 2

function [emd_mat,eb_loc] = build_emd(N),

%total number of edges
te = 2*N - 1;

%find length of message communicated over each edge
de = floor(log2(1:te));
mle = N ./ (2.^de);

%eb_loc finds the storage location in the path for DWM and UWM 
% transmitted over edges (column 1 for DWM, column 2 for UWM)
eb_loc = zeros(2,te);

%initially a DWM is passed through edge 1
%This intial move is not stored in emd_mat
CMD = 1; CE = 1;

%the initial path length is the length of channel beliefs
%which is equal to the codeword length (N)
cpl = N;

%save the start location of the intial message (which is location 1 in the
%path)
eb_loc(1,1) = 1;

emd_mat = [];


while (1), 

    if CMD == 2 && CE == 1, %case 5: termination of process
        break;
    
    elseif CMD == 1 && CE >= N, %case 1: DWM over a LIFE
        
        NMD = 2; NE = CE; 

    elseif CMD == 1 && CE < N, %case 2: DWM on NLIF

        NMD = 1; NE = 2*CE;

    elseif CMD == 2 && mod(CE,2)  == 0; %case 3: UWM on LCE

        NMD = 1; NE = CE+1;

    elseif CMD == 2 && mod(CE,2) == 1, %case 4: UMW on RCE

        NMD = 2; NE = floor(CE/2);

    else
        disp('Unexpected case in function build_emd');
    end

    %%%% update matrices and path length
    CE = NE; CMD = NMD;

    emd_mat = [emd_mat, [CE; CMD]];

    eb_loc(CMD,CE) = cpl+1;

    cpl = cpl + mle(CE);

end %end of while loop




return