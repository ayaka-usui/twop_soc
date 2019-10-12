function f = C2u(u,ksoc)
    mat = zeros(4,4);
    mat(4,2) = -1i*ksoc*sqrt(u);
    f = mat;
end