function f = matF3_ver2(N,M,fuv,g,g12) % N<=M
    matF2 = zeros(4*(N+1),4*(N+1));
    for nn = 1:N+1
        matF0 = zeros(N+1,N+1);
        matF0(nn,nn) = 1;
        matF1 = kron(matF0,matF(fuv(N+1-nn+1,M+1-nn+1),g,g12));        
        matF2 = matF2 + matF1;
    end
    matF2 = sparse(matF2);
    f = horzcat(matF2,sparse(4*(N+1),4*(M-N)));
end