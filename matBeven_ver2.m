function f = matBeven_ver2(N,ksoc)
    matB2even = zeros(2*N,2*N);
    for nn = 1:N/2
        matB0even = zeros(N/2,N/2);
        matB0even(nn,nn) = 1;
        matB1even = kron(matB0even,Bn(2*nn,ksoc));
        matB2even = matB2even + matB1even;
    end
    matB2even = sparse(matB2even);
    matB2even = horzcat(sparse(2*N,4),matB2even);
    f = matB2even;
end