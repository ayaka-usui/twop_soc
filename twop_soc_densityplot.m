
% Load the data produced "data_eigEVgs_Omegaj_k' num2str(ksoc) '_g' num2str(g) '_g12' num2str(g12) '_Nsize' num2str(Nsize) '.mat"
% This file makes Fig. 4 and 7
% From Line 179, a section for plot starts. Run a section where you want to have a plot

%% Parameters
N0 = 99;
Nsize = 4*(N0+1)*(N0+2);

g = 0; %0.4; %0
g12 = 10; %1; %0; %10
ksoc = 5; %1.5;

%% Eigenstates of harmonic oscillator
dx = 0.1; %0.01; %0.05;
Xmax = 20; %10; %5;
x = (-Xmax:dx:Xmax);
Nx = length(x);

pxmax = pi*Nx/(2*Xmax);
dpx = 2*pxmax/Nx;
px = ((1:Nx)*dpx-pxmax);

psiH = zeros(Nx,Nsize);
for nn = 0:N0
    psiH(:,nn+1) = hermiteH(nn,x(:)).*exp(-x(:).^2/2);
    psiH(:,nn+1) = psiH(:,nn+1)/sqrt(sum(abs(psiH(:,nn+1)).^2)*dx);
end

% save(['data_psiH_Nsize' num2str(Nx) 'dx' num2str(dx) 'Xmax' num2str(Xmax) '.mat'], ...
%       'x','psiH')

%% Index for matrices of an
inddown = (1:4:Nsize-3);
indS = (2:4:Nsize-2);
indup = (3:4:Nsize-1);
indA = (4:4:Nsize);
indsize = length(inddown);

indmat=zeros(indsize,2);
jj = 0;
for nn = 0:Nsize-1
    for mm = 0:nn
        if jj >= indsize
           break
        end
        jj = jj + 1;
        indmat(jj,:) = [2*mm; 2*(nn-mm)];
    end
    for mm = 0:nn
        if jj >= indsize
           break
        end
        jj = jj + 1;
        indmat(jj,:) = [2*mm+1; 2*(nn-mm)];
    end
end
indmat = indmat + ones(indsize,2);

%% Density 
Nomega = length(Omegaj);

phidownj = zeros(Nx,Nx,Nomega);
% phiSj = zeros(Nx,Nx,Nomega);
phiupj = zeros(Nx,Nx,Nomega);
% phiAj = zeros(Nx,Nx,Nomega);
phidownupj = zeros(Nx,Nx,Nomega);
phiupdownj = zeros(Nx,Nx,Nomega);
densityj = zeros(Nx,Nx,Nomega);

phidownjF = zeros(Nx,Nx,Nomega);
% phiSjF = zeros(Nx,Nx,Nomega);
phiupjF = zeros(Nx,Nx,Nomega);
% phiAjF = zeros(Nx,Nx,Nomega);
phidownupjF = zeros(Nx,Nx,Nomega);
phiupdownjF = zeros(Nx,Nx,Nomega);
densityjF = zeros(Nx,Nx,Nomega);

for jj = 60:100 %1:Nomega 
    %jj=60:100 means Omegaj=30~50
    %jj=1:Nomega takes long time, and so you may want to pick up only a range you need
    
    an = eigVgs(:,jj);
    
    andown = an(inddown);
    anS = an(indS);
    anup = an(indup);
    anA = an(indA);

    % mat_an
    matandown = zeros(1,1);
    for nn = 1:indsize
        matandown(indmat(nn,1),indmat(nn,2)) = andown(nn);
    end
    matanS = zeros(1,1);
    for nn = 1:indsize
        matanS(indmat(nn,1),indmat(nn,2)) = anS(nn);
    end
    matanup = zeros(1,1);
    for nn = 1:indsize
        matanup(indmat(nn,1),indmat(nn,2)) = anup(nn);
    end
    matanA = zeros(1,1);
    for nn = 1:indsize
        matanA(indmat(nn,1),indmat(nn,2)+1) = anA(nn); % Note that anA = anA(int,odd) while anS = anS(int,even)
    end
    clear andown anS anup anA
    
    % Density distribution
    nn_length0 = size(matandown,1)-1;
    mm_length0 = floor((size(matandown,2)-1)/2);
    
    phidown = zeros(Nx,Nx);
    for nn = 0:nn_length0
        for mm = 0:mm_length0
            phidown = phidown + matandown(nn+1,2*mm+1)*psiH(:,nn+1)*psiH(:,2*mm+1).';
        end
    end
    
    phiS = zeros(Nx,Nx);
    for nn = 0:nn_length0
        for mm = 0:mm_length0
            phiS = phiS + matanS(nn+1,2*mm+1)*psiH(:,nn+1)*psiH(:,2*mm+1).';
        end
    end
    
    phiup = zeros(Nx,Nx);
    for nn = 0:nn_length0
        for mm = 0:mm_length0
            phiup = phiup + matanup(nn+1,2*mm+1)*psiH(:,nn+1)*psiH(:,2*mm+1).';
        end
    end
    
    phiA = zeros(Nx,Nx);
    for nn = 0:nn_length0
        for mm = 0:mm_length0
            phiA = phiA + matanA(nn+1,2*mm+1+1)*psiH(:,nn+1)*psiH(:,2*mm+1+1).';
        end
    end

    phidownup = (phiS - phiA)/sqrt(2);
    phiupdown = (phiS + phiA)/sqrt(2);
    
    % Momentum distribution (FFT)
    phidownjF(:,:,jj) = fftshift(fft2(phidown));
%     phiSjF(:,:,jj) = fftshift(fft2(phiS));
    phiupjF(:,:,jj) = fftshift(fft2(phiup));
%     phiAjF(:,:,jj) = fftshift(fft2(phiA));
    phidownupjF(:,:,jj) = fftshift(fft2(phidownup));
    phiupdownjF(:,:,jj) = fftshift(fft2(phiupdown));
    
    % Normalisation
    % for density distribution
    norm = sqrt(sum(sum(abs(phidown).^2 + abs(phidownup).^2 + abs(phiup).^2 + abs(phiupdown).^2))*dx*dx);
    phidownj(:,:,jj) = phidown/norm;
%     phiSj(:,:,jj) = phiS/norm;
    phiupj(:,:,jj) = phiup/norm;
%     phiAj(:,:,jj) = phiA/norm;
    phidownupj(:,:,jj) = phidownup/norm;
    phiupdownj(:,:,jj) = phiupdown/norm;
    
    % for momentum distribution
    normF = sqrt(sum(sum(abs(phidownjF(:,:,jj)).^2 + abs(phidownupjF(:,:,jj)).^2 + abs(phiupjF(:,:,jj)).^2 + abs(phiupdownjF(:,:,jj)).^2))*dpx*dpx);
    phidownjF(:,:,jj) = phidownjF(:,:,jj)/normF;
    phiupjF(:,:,jj) = phiupjF(:,:,jj)/normF;
    phidownupjF(:,:,jj) = phidownupjF(:,:,jj)/normF;
    phiupdownjF(:,:,jj) = phiupdownjF(:,:,jj)/normF;
    
end

% save(['data_densitymomentum_correlation_k' num2str(ksoc) '_g' num2str(g) 'g12' num2str(g12) '_Nsize' num2str(Nx) 'dx' num2str(dx) 'Xmax' num2str(Xmax) '.mat'], ...
%       'x','phidownj','phidownupj','phiupdownj','phiupj', ...
%       'px','phidownjF','phidownupjF','phiupdownjF','phiupjF', ...
%       'densityj','densityjF')

%%
% Run a section where you want to have a plot

%% Momentum distribution of each pseudo-spin state, downdown, downup, updown, upup

figure
for jj = 60:100 %1:Nomega
    
    subplot(2,2,1)
    pcolor(px,px,abs(phidownjF(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    title(['\Omega=' num2str(Omegaj(jj))])
    xlim([-10 10])
    ylim([-10 10])
    % caxis([0 0.22])
    % caxis([0 0.3206])
    % caxis([0 0.1445])
    % caxis([0 0.3393])
    daspect([1 1 1])
    
    subplot(2,2,2)
    pcolor(px,px,abs(phidownupjF(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    xlim([-10 10])
    ylim([-10 10])
    % caxis([0 0.22])
    % caxis([0 0.3206])
    % caxis([0 0.1445])
    % caxis([0 0.3393])
    daspect([1 1 1])
    
    subplot(2,2,3)
    pcolor(px,px,abs(phiupdownjF(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    xlim([-10 10])
    ylim([-10 10])
    % caxis([0 0.22])
    % caxis([0 0.3206])
    % caxis([0 0.1445])
    % caxis([0 0.3393])
    daspect([1 1 1])

    subplot(2,2,4)
    pcolor(px,px,abs(phiupjF(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    xlim([-10 10])
    ylim([-10 10])
    % caxis([0 0.22])
    % caxis([0 0.3206])
    % caxis([0 0.1445])
    % caxis([0 0.3393])
    daspect([1 1 1])
    
    pause %(0.01)

end

%% Overall momentum distribution, Fig. 4(a-d) and 7(a-d)

figure
for jj = 60:100 %1:Nomega
    
    pcolor(px,px,abs(phidownjF(:,:,jj))+abs(phidownupjF(:,:,jj))+abs(phiupdownjF(:,:,jj))+abs(phiupjF(:,:,jj)))
    set(gca,'FontSize',32,'FontName','Times New Roman')
    shading flat
    colorbar
    title(['\Omega=' num2str(Omegaj(jj))])
    xlim([-10 10])
    ylim([-10 10])
    daspect([1 1 1])
    
    pause
    
end

%% Density distribution of each pseudo-spin state, downdown, downup, updown, upup

figure
for jj = 60:100 %1:Nomega

    subplot(2,2,1)
    pcolor(x,x,abs(phidownj(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    title(['\Omega=' num2str(Omegaj(jj))])
    xlim([-3 3])
    ylim([-3 3])
    % caxis([0 0.57])
    daspect([1 1 1])

    subplot(2,2,2)
    pcolor(x,x,abs(phidownupj(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    xlim([-3 3])
    ylim([-3 3])
    % caxis([0 0.57])
    daspect([1 1 1])

    subplot(2,2,3)
    pcolor(x,x,abs(phiupdownj(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    xlim([-3 3])
    ylim([-3 3])
    % caxis([0 0.57])
    daspect([1 1 1])

    subplot(2,2,4)
    pcolor(x,x,abs(phiupj(:,:,jj)))
    set(gca,'FontSize',20,'FontName','Times New Roman')
    shading flat
    colorbar
    xlim([-3 3])
    ylim([-3 3])
    % caxis([0 0.57])
    daspect([1 1 1])

    pause %(0.01)

end

%% Overall density distribution, Fig. 4(e-h) and 7(e-h)

figure
for jj = 60:100 %1:Nomega
    
    pcolor(px,px,abs(phidownj(:,:,jj))+abs(phidownupj(:,:,jj))+abs(phiupdownj(:,:,jj))+abs(phiupj(:,:,jj)))
    set(gca,'FontSize',32,'FontName','Times New Roman')
    shading flat
    colorbar
    title(['\Omega=' num2str(Omegaj(jj))])
    xlim([-5 5])
    ylim([-5 5])
    daspect([1 1 1])
    
    pause
    
end
