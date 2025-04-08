% x = polar_encode(m,N)
%x is the codeword, m is the input bit sequence, N is the codeword length
function x = polar_encode(m,N,G,Rseq),

n = floor(log2(N));

K = length(m);

if (K > N || n ~= log2(N)),
    disp('wrong inputs to function polar_encode');
    return;
end

% %Load the 5G standard reliability sequence for 1024
% rel_seq1024;
% 
% %Extract reliability sequence for N
% [ir,jr] = find(Rseq1024 <= N);
% 
% Rseq = Rseq1024(jr);
% 
% %Build polar transform matrix for n
% G = build_G(n);

u = zeros(1,N);

u(Rseq(N-K+1:N)) = m;

x = rem(u*G,2);

return