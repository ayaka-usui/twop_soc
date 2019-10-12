function f = Bn(n,ksoc)
    mat = zeros(4,4);
    mat(1,1) = 1i*ksoc*sqrt(n);
    mat(3,3) =-1i*ksoc*sqrt(n);
    f = mat;
end