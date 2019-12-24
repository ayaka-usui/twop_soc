function f = matAeven_ver2(N,ksoc,Omega)
    matA2even = zeros(2*N+4,2*N+4);
    for nn = 1:N/2+1
        matA0even = zeros(N/2+1,N/2+1);
        matA0even(nn,nn) = 1;
        matA1even = kron(matA0even,An2u(2*(nn-1),N-2*(nn-1),ksoc,Omega));        
        matA2even = matA2even + matA1even;
    end
    f = sparse(matA2even);
end
