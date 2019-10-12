clear all
tic

%% Parameters
g = 0;
g12 = 1;
ksoc = 5;

N0 = 99;
Nsize = 4*(N0+1)*(N0+2);
matH = spalloc(Nsize,Nsize,round(Nsize*Nsize*0.1*0.05*0.4));
Omegaj = (0.5:0.5:60);
Nomega = length(Omegaj);
Nspec = 4*2*3;

%% Noninteracting part independent of Omega
for ss = 0:N0
    % matB
    matH(4*ss*(ss+1)+1:4*(ss+1)^2,4*(ss+1)^2+1:4*(ss+1)*(ss+2)) = matBodd_ver2(2*ss+1,ksoc);
    matH(4*(ss+1)^2+1:4*(ss+1)*(ss+2),4*ss*(ss+1)+1:4*(ss+1)^2) = matBodd_ver2(2*ss+1,ksoc)';
    if ss == 0
       continue % ss>=1
    end
    matH(4*ss^2+1:4*ss*(ss+1),4*ss*(ss+1)+1:4*(ss+1)^2) = matBeven_ver2(2*ss,ksoc);
    matH(4*ss*(ss+1)+1:4*(ss+1)^2,4*ss^2+1:4*ss*(ss+1)) = matBeven_ver2(2*ss,ksoc)';
    
    % matC
    matH(4*ss^2+1-4*ss:4*ss*(ss+1)-4*ss,4*ss*(ss+1)+1:4*(ss+1)^2) = matCeven_ver2(2*ss,ksoc);
    matH(4*ss*(ss+1)+1:4*(ss+1)^2,4*ss^2+1-4*ss:4*ss*(ss+1)-4*ss) = matCeven_ver2(2*ss,ksoc)';
    matH(4*ss^2+1:4*ss*(ss+1),4*ss*(ss+1)+1+4*(ss+1):4*(ss+1)^2+4*(ss+1)) = matCeven_ver2(2*ss,ksoc);
    matH(4*ss*(ss+1)+1+4*(ss+1):4*(ss+1)^2+4*(ss+1),4*ss^2+1:4*ss*(ss+1)) = matCeven_ver2(2*ss,ksoc)';
end
matH = sparse(matH);

%% Interaction matF
% fuv = f(u+1,v+1)
f00 = 1/sqrt(2*pi);
indmax = N0;
fuv = zeros(indmax+1,indmax+1);
for uu = 0:indmax
    for vv = uu:indmax
        if vv==0
           fuv(uu+1,vv+1) = f00;
        elseif uu > vv-1
           fuv(uu+1,vv+1) = -sqrt((2*vv-1)/2/vv)*fuv(vv,uu+1);
        else
           fuv(uu+1,vv+1) = -sqrt((2*vv-1)/2/vv)*fuv(uu+1,vv);
        end 
    end
    fuv(:,uu+1) = fuv(uu+1,:).';
end

%%
g12array = (0:0.25:10);
for jg = 1:length(g12array)

    g12 = g12array(jg);
    
    % matHint
    % off-diagonal terms
    matHint = spalloc(Nsize,Nsize,round(Nsize*Nsize*0.1*0.05*0.4)); %zeros(Nsize,Nsize);
    for jj = 1:N0
        for ss = jj:N0
            matHint(4*(ss-jj)*(ss-jj+1)+1:4*(ss-jj+1)^2,4*ss*(ss+1)+1:4*(ss+1)^2) = matF3_ver2(ss-jj,ss,fuv(1:ss-jj+1,1:ss+1),g,g12);
            matHint(4*(ss-jj+1)^2+1:4*(ss-jj+1)*(ss-jj+2),4*(ss+1)^2+1:4*(ss+1)*(ss+2)) = matF3_ver2(ss-jj,ss,fuv(1:ss-jj+1,1:ss+1),g,g12);
        end
    end
    matHint = matHint + matHint';

    % diagonal terms
    for ss = 0:N0
        matHint(4*ss*(ss+1)+1:4*(ss+1)^2,4*ss*(ss+1)+1:4*(ss+1)^2) = matF3_ver2(ss,ss,fuv(1:ss+1,1:ss+1),g,g12);
        matHint(4*(ss+1)^2+1:4*(ss+1)*(ss+2),4*(ss+1)^2+1:4*(ss+1)*(ss+2)) = matF3_ver2(ss,ss,fuv(1:ss+1,1:ss+1),g,g12);
    end

    matHint = sparse(matHint);

    % Diagonalisation
    matH1 = matH + matHint;

    % dependent part of Omega
    Espec = zeros(Nspec,Nomega);
    eigVgs = zeros(Nsize,Nomega);
    for jj = 1:Nomega
        matHA = spalloc(Nsize,Nsize,round(Nsize*Nsize*0.1*0.05*0.4)); %zeros(Nsize,Nsize);
        Omega = Omegaj(jj);
        for ss = 0:N0
            % matA
            matHA(4*ss*(ss+1)+1:4*(ss+1)^2,4*ss*(ss+1)+1:4*(ss+1)^2) = matAeven_ver2(2*ss,ksoc,Omega);
            matHA(4*(ss+1)^2+1:4*(ss+1)*(ss+2),4*(ss+1)^2+1:4*(ss+1)*(ss+2)) = matAodd_ver2(2*ss+1,ksoc,Omega);
        end
        matHA = sparse(matHA);
        matHT = matH1 + matHA;
    
        [eigV,eigE] = eigs(matHT,Nspec,'sr');
        [eigE,Iorder] = sort(diag(eigE));
        Espec(:,jj) = eigE(:);
        eigV = eigV(:,Iorder);
        eigVgs(:,jj) = eigV(:,1);
        
    end

    % Save data
    if 1
    cd data
    save(['data_eigEVgs_Omegaj_k' num2str(ksoc) '_g' num2str(g) '_g12' num2str(g12) '_Nsize' num2str(Nsize) '.mat'], ...
          'Omegaj','Nomega','Espec','eigVgs','Nsize')
    cd ..
    end

end

%%
toc















