function f = An2u(n,u,ksoc,Omega)
    mat = zeros(4,4);
    
    mat(1,2) = Omega/sqrt(2);
    mat(2,3) = Omega/sqrt(2);
    mat(2,4) = -1i*ksoc*sqrt(u+1);
    mat = mat + mat';
    
    for jj = 1:4
        mat(jj,jj) = n + u + 1;
        %
        mat(jj,jj) = mat(jj,jj) + ksoc^2;
    end
    mat(4,4) = mat(4,4) + 1;
    
    f = mat;
end