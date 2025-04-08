This simulation package includes MATLAB codes developed for the following work:

[1] J. Haghighat and T. M. Duman, “A Practical Indexing Scheme for Noisy Shuffling Channels Using Cosets of Polar Codes,” IEEE Transactions on Communications, doi: 10.1109/TCOMM.2024.3506944.

Main MATLAB Codes:

RS_Polar.m
Simulates a concatenated Reed–Solomon (RS) and Polar coding/decoding scheme over a noisy shuffling channel.

RS_LDPC.m
Simulates a concatenated Reed–Solomon (RS) and LDPC coding/decoding scheme over a noisy shuffling channel.

polar_list_dec_BSC_v3.m and build_emd.m
Implement the list decoder for polar codes based on the following reference:
I. Tal and A. Vardy, “List decoding of polar codes,” IEEE Transactions on Information Theory, vol. 61, no. 5, pp. 2213–2226, May 2015.

Supporting Functions:

polar_encode.m
Implements the polar encoding process.

build_G.m
Generates the polar transform matrix.

rel_seq1024.m
Contains the 5G standard channel-reliability sequence used in polar coding, as specified in:
V. Bioglio, C. Condo, and I. Land, “Design of polar codes in 5G New Radio,” IEEE Communications Surveys & Tutorials, vol. 23, no. 1, pp. 29–40, Q1 2021.

sort_metric_opt.m and sort_metric.m
Implement optimal and suboptimal sorting algorithms for data received through the shuffling channel. These are detailed in [1].

build_incomplete_Gallager.m
Generates the parity-check and generator matrices of an LDPC code using Gallager’s construction.

ldpc_dec_sp.m
Implements the sum-product algorithm for decoding noisy LDPC codewords.

For further inquiries, please contact Javad Haghighat at:
javad.haghighat@tedu.edu.tr

