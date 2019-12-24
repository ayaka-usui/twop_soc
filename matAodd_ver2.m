function f = matAodd_ver2(N,ksoc,Omega)
    matA2odd = zeros(2*N+2,2*N+2);
    for nn = 1:(N+1)/2
        matA0odd = zeros((N+1)/2,(N+1)/2);
        matA0odd(nn,nn) = 1;
        matA1odd = kron(matA0odd,An2u(2*nn-1,N-1-2*(nn-1),ksoc,Omega));        
        matA2odd = matA2odd + matA1odd;
    end
    f = sparse(matA2odd);
end
