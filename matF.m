function f = matF(fuv,g,g12)
    mat = zeros(4,4);
    mat(1,1) = g*fuv;
    mat(2,2) = g12*fuv;
    mat(3,3) = g*fuv;
%     matF(4,4) = 0; 
    f = mat;
end