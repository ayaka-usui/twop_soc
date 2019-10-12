function f = matBodd_ver2(N,ksoc)
    matB2odd = zeros(2*N+2,2*N+2);
    for nn = 1:(N+1)/2
        matB0odd = zeros((N+1)/2,(N+1)/2);
        matB0odd(nn,nn) = 1;
        matB1odd = kron(matB0odd,Bn(2*nn-1,ksoc));
        matB2odd = matB2odd + matB1odd;
    end
    f = sparse(matB2odd);
end