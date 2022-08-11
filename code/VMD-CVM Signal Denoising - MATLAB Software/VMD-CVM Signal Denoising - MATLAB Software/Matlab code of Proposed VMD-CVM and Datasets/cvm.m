function test=cvm(z,win_len)
N=win_len;
I=1:N;
test = sum((z-(2*I-1)/2/N).^2) + 1/12/N;
end