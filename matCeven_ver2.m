function f = matCeven_ver2(N,ksoc)
    matC2even = zeros(2*N,2*N);
    for nn = 1:N/2
        matC0even = zeros(N/2,N/2);
        matC0even(N/2+1-nn,N/2+1-nn) = 1;
        matC1even = kron(matC0even,C2u(2*nn,ksoc));
        matC2even = matC2even + matC1even;
    end
    matC2even = sparse(matC2even);
    matC2even = horzcat(matC2even,sparse(2*N,4));
    f = matC2even;
end