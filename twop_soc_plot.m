
% Load the data produced "data_eigEVgs_Omegaj_k' num2str(ksoc) '_g' num2str(g) '_g12' num2str(g12) '_Nsize' num2str(Nsize) '.mat"
% Run a section where you want to have a plot
% This file makes all the figures except for Fig. 4 and 7

%% Energy spectrum, Fig. 2 and 6(a,b)
figure
for nn = 2:Nspec
    plot(Omegaj(:),Espec(nn,:)-Espec(1,:),'-','LineWidth',3)
    hold on
end
ylim([0 1.1])
set(gca,'FontSize',32,'FontName','Times New Roman')

%% Population, Fig. 3(a,b) and 6(c)
antot = zeros(length(Omegaj),4);
for jj = 1:length(Omegaj)
    an = eigVgs(:,jj);
    for nn = 0:Nsize/4-1
        antot(jj,1) = antot(jj,1) + abs(an(4*nn+1)).^2;
        antot(jj,2) = antot(jj,2) + abs(an(4*nn+2)).^2;
        antot(jj,3) = antot(jj,3) + abs(an(4*nn+3)).^2;
        antot(jj,4) = antot(jj,4) + abs(an(4*nn+4)).^2;
    end
end

% Population of downdown, S, upup, and A
figure
plot(Omegaj,antot(:,1),'-o')
hold on
plot(Omegaj,antot(:,2),'-s')
plot(Omegaj,antot(:,3),'-x')
plot(Omegaj,antot(:,4),'-*')
legend('downdown','S','upup','A')

% Population of downdown, downup, updown, upup
figure
plot(Omegaj,antot(:,1),'LineStyle','-','Color',[0 0 0],'LineWidth',3)
hold on
plot(Omegaj,(antot(:,2)+antot(:,4))/2,'LineStyle','-','Color',[0 0 1],'LineWidth',3)
plot(Omegaj,(antot(:,2)+antot(:,4))/2,'LineStyle',':','Color',[0 1 1],'LineWidth',3)
plot(Omegaj,antot(:,3),'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',3)
legend('\downarrow\downarrow','\downarrow\uparrow','\uparrow\downarrow','\uparrow\uparrow')
set(gca,'FontSize',32,'FontName','Times New Roman')

%% Interaction energy, Fig. 3(a,b)
g = 0;
g12 = 1;
Eintj = zeros(Nomega,1);

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

for jj=1:Nomega
    an = eigVgs(:,jj);
    Eintj(jj) = an'*(matHint*an);
end

figure
plot(Omegaj,Eintj,'-k','LineWidth',3)
set(gca,'FontSize',32,'FontName','Times New Roman')

%% Population difference, Fig. 3(c)
g12array = (0:0.25:10);
antot_g12 = zeros(length(Omegaj),4,length(g12array));
for jg = 1:length(g12array)
    g12 = g12array(jg);
    
    % Load data
    load(['/Volumes/work/BuschU/ayakausui/simon_diag/g12varry/k5/data/data_eigEVgs_Omegaj_k5_g0_g12' num2str(g12) '_Nsize40400.mat'])
    
    antot = zeros(length(Omegaj),4);
    for jj = 1:length(Omegaj)
        an = eigVgs(:,jj);
        for nn = 0:Nsize/4-1
            antot(jj,1) = antot(jj,1) + abs(an(4*nn+1)).^2;
            antot(jj,2) = antot(jj,2) + abs(an(4*nn+2)).^2;
            antot(jj,3) = antot(jj,3) + abs(an(4*nn+3)).^2;
            antot(jj,4) = antot(jj,4) + abs(an(4*nn+4)).^2;
        end
    end
    
    antot_g12(:,:,jg) = antot(:,:);
    
end

clear diff
for jw = 1:length(Omegaj)
    for jg = 1:length(g12array)
        diff(jw,jg) = (antot_g12(jw,2,jg)+antot_g12(jw,4,jg))/2-antot_g12(jw,1,jg);
        if diff(jw,jg) <= 0 || jg == 1
           diff(jw,jg) = NaN;
        end
    end
end

figure
pcolor(Omegaj,g12array,diff')
xlim([35 60])
set(gca,'FontSize',32,'FontName','Times New Roman')
shading flat
colorbar

%% Population difference, Fig. 6(d)
garray = (0:0.02:2);
antot_g = zeros(length(Omegaj),4,length(garray));
for jg = 1:length(garray)
    g = garray(jg);
    
    % Load data
%     load(['/Volumes/work/BuschU/ayakausui/simon_diag/gvarry/g_01_3/data/data_eigEVgs_Omegaj_k5_g' num2str(g) '_g120_Nsize40400.mat'])
%     load(['/Volumes/work/BuschU/ayakausui/simon_diag/g00_025_3_up/data/data_eigEVgs_Omegaj_k5_g' num2str(g) '_g120_Nsize40400.mat'])
    load(['/Volumes/work/BuschU/ayakausui/simon_diag/gvarry/g0_0.02_2/data/data_eigEVgs_Omegaj_k5_g' num2str(g) '_g120_Nsize40400.mat'])
    
    antot = zeros(length(Omegaj),4);
    for jj = 1:length(Omegaj)
        an = eigVgs(:,jj);
        for nn = 0:Nsize/4-1
            antot(jj,1) = antot(jj,1) + abs(an(4*nn+1)).^2;
            antot(jj,2) = antot(jj,2) + abs(an(4*nn+2)).^2;
            antot(jj,3) = antot(jj,3) + abs(an(4*nn+3)).^2;
            antot(jj,4) = antot(jj,4) + abs(an(4*nn+4)).^2;
        end
    end
    
    antot_g(:,:,jg) = antot(:,:);
    
end

diff = zeros(length(Omegaj),length(garray));
for jw = 1:length(Omegaj)
    for jg = 1:length(garray)
        diff(jw,jg) = antot_g(jw,1,jg) - (antot_g(jw,2,jg)+antot_g(jw,4,jg))/2;
        if diff(jw,jg) <= 0 || jg == 1
           diff(jw,jg) = NaN;
        end
    end
end

figure
pcolor(Omegaj,garray,diff')
set(gca,'FontSize',32,'FontName','Times New Roman')
shading flat
colorbar

%% Coeffieient an of Eq. (6) (This is needed for the below plots)
inddown = (1:4:Nsize-3);
indS = (2:4:Nsize-2);
indup = (3:4:Nsize-1);
indA = (4:4:Nsize);
indsize = length(inddown);

% index for matrices of an
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

%% Evaluation of Entanglement, Fig. 5(c) and 8(c)
% Concurrence and von Neumann entropy as functions of Omega

rhoj = zeros(4,4,Nomega);
Concj = zeros(Nomega,1);
Concj1 = zeros(Nomega,1);
lamda1j = zeros(Nomega,4);
vec1j = zeros(Nomega,4,4); % should be vec1j = zeros(4,4,Nomega); though
vNEj = zeros(Nomega,1);
rhoconj = zeros(4,4,Nomega);
lamdaj = zeros(Nomega,4);
linentroj = zeros(Nomega,1);
rhoPTj = zeros(4,4,Nomega);
lamdaPTj = zeros(Nomega,1);

for jj = 1:Nomega
    % an
    an = eigVgs(:,jj);
    andown = an(inddown);
    anS = an(indS);
    anup = an(indup);
    anA = an(indA);
    
    rho = zeros(4,4);
    
    rho(1,1) = sum(abs(andown(:)).^2);
    rho(1,2) = sum(conj(andown(:)).*anS(:))/sqrt(2);
    rho(1,3) = sum(conj(andown(:)).*anS(:))/sqrt(2);
    rho(1,4) = sum(conj(andown(:)).*anup(:));

    rho(2,2) = sum(abs(anS(:)).^2 + abs(anA(:)).^2)/2;
    rho(2,3) = sum(abs(anS(:)).^2 - abs(anA(:)).^2)/2;
    rho(2,4) = sum(conj(anS(:)).*anup(:))/sqrt(2);

    rho(3,3) = sum(abs(anS(:)).^2 + abs(anA(:)).^2)/2;
    rho(3,4) = sum(conj(anS(:)).*anup(:))/sqrt(2);

    rho(4,4) = sum(abs(anup(:)).^2);

    rho = rho + (rho - diag(diag(rho)))';
    rho = rho/trace(rho);
    
    rho = real(rho);
    rhoj(:,:,jj) = rho;
            
    % Concurrence
    sigmay = [0 -1i; 1i 0];
    sigmay2 = kron(sigmay,sigmay);
%     lamda = eig(rho*sigmay2*conj(rho)*sigmay2);
    lamda2 = eig(sqrtm(sqrtm(rho)*sigmay2*conj(rho)*sigmay2*sqrtm(rho)));
%     lamda0 = lamda;
%     lamda = sort(lamda,'descend'); 
    lamda = sort(lamda2,'descend');
    Conc = lamda(1)-lamda(2)-lamda(3)-lamda(4); %sqrt(lamda(1))-sqrt(lamda(2))-sqrt(lamda(3))-sqrt(lamda(4)); %lamda(1)-lamda(2)-lamda(3)-lamda(4);
    Concj(jj) = Conc;
    if Conc <= 0
       Conc = 0;
    end
    Concj1(jj) = Conc;
    
    rhoconj(:,:,jj) = rho*sigmay2*conj(rho)*sigmay2;
    lamdaj(jj,:) = lamda; %lamda0;
    
    % von Neumman entropy
    [vec1,lamda1] = eig(rho);
    lamda1 = diag(lamda1);
    [lamda1,ind] = sort(lamda1);
    vec1j(jj,:,:) = vec1(:,ind);
    vNE = 0;
    for nn = 1:4
        if lamda1(nn) ~= 0
           vNE = vNE - sum(lamda1(nn)*log2(lamda1(nn)));
        end
    end
    lamda1j(jj,:) = lamda1(:);
    vNEj(jj) = vNE;
   
    % Linear entropy
    linentroz = 1 - trace(rho^2);
    linentroj(jj) = linentroz;
    
    % Partial transport to detect entanglement
    rho1 = rho(1:2,1:2).';
    rho2 = rho(1:2,3:4).';
    rho3 = rho(3:4,1:2).';
    rho4 = rho(3:4,3:4).';
    
    rhoPT1 = [rho1 rho2];
    rhoPT2 = [rho3 rho4];
    rhoPT = [rhoPT1; rhoPT2];
    lamdaPT0 = eig(rhoPT);
    indPT = find(lamdaPT0<0);
    lamdaPT1 = -sum(lamdaPT0(indPT));
    
    rhoPTj(:,:,jj) = rhoPT;
    lamdaPTj(jj) = lamdaPT1;
    
end

% von Neumann entropy
figure
plot(Omegaj(:),lamda1j(:,1),':','LineWidth',3)
hold on
plot(Omegaj(:),lamda1j(:,2),':','LineWidth',3)
plot(Omegaj(:),lamda1j(:,3),':','LineWidth',3)
plot(Omegaj(:),lamda1j(:,4),':','LineWidth',3)
plot(Omegaj,vNEj,'-k','LineWidth',3)
set(gca,'FontSize',32,'FontName','Times New Roman')

%% Linear Entropy
figure
plot(Omegaj,linentroj,'-o')

%% Partial transport to detect entanglement
figure
plot(Omegaj,lamdaPTj,'-o')

%% Concurrence
figure
plot(Omegaj,real(Concj1),'-k','LineWidth',3) %real(Concj1)
% ylim([0 4*10^(-3)])
set(gca,'FontSize',32,'FontName','Times New Roman')

% Lamda shown in Eq. (10)
figure
plot(Omegaj(:),lamdaj(:,1),'-','LineWidth',3)
hold on
plot(Omegaj(:),lamdaj(:,2),'-','LineWidth',3)
plot(Omegaj(:),lamdaj(:,3),'-','LineWidth',3)
plot(Omegaj(:),lamdaj(:,4),'-','LineWidth',3)
set(gca,'FontSize',32,'FontName','Times New Roman')

%% Concurrence and von Neumann entropy as functions of ksoc and Omega, Fig. 5(a,b) and 8(a,b)
ksocarray = (5:-0.1:0); %(5:-0.25:0);
Concj_ksoc = zeros(Nomega,length(ksocarray));
vNEj_ksoc = zeros(Nomega,length(ksocarray));

for jk = 1:length(ksocarray)
    ksoc = ksocarray(jk);
    
    % Load data
    load(['/Volumes/work/BuschU/ayakausui/simon_diag/varryk/g04/data/data_eigEVgs_Omegaj_k' num2str(ksoc) '_g0.4_g120_Nsize40400.mat'])

%     Concj = zeros(Nomega,1);
    Concj1 = zeros(Nomega,1);

    for jj = 1:Nomega
        % an
        an = eigVgs(:,jj);
        andown = an(inddown);
        anS = an(indS);
        anup = an(indup);
        anA = an(indA);
        
        rho = zeros(4,4);
    
        rho(1,1) = sum(abs(andown(:)).^2);
        rho(1,2) = sum(conj(andown(:)).*anS(:))/sqrt(2);
        rho(1,3) = sum(conj(andown(:)).*anS(:))/sqrt(2);
        rho(1,4) = sum(conj(andown(:)).*anup(:));
        
        rho(2,2) = sum(abs(anS(:)).^2 + abs(anA(:)).^2)/2;
        rho(2,3) = sum(abs(anS(:)).^2 - abs(anA(:)).^2)/2;
        rho(2,4) = sum(conj(anS(:)).*anup(:))/sqrt(2);
        
        rho(3,3) = sum(abs(anS(:)).^2 + abs(anA(:)).^2)/2;
        rho(3,4) = sum(conj(anS(:)).*anup(:))/sqrt(2);
        
        rho(4,4) = sum(abs(anup(:)).^2);
        
        rho = rho + (rho - diag(diag(rho)))';
        rho = rho/trace(rho);
            
        % Concurrence
        sigmay = [0 -1i; 1i 0];
        sigmay2 = kron(sigmay,sigmay);
        lamda = eig(rho*sigmay2*conj(rho)*sigmay2);
        lamda = sort(lamda,'descend');
        Conc = sqrt(lamda(1))-sqrt(lamda(2))-sqrt(lamda(3))-sqrt(lamda(4)); %lamda(1)-lamda(2)-lamda(3)-lamda(4);
%         Concj(jj) = Conc;
%         if Conc <= 0
%            Conc = 0;
%         end
        Concj1(jj) = Conc;
        
        % von Neumman entropy
        [~,lamda1] = eig(rho); %[vec1,lamda1] = eig(rho);
        lamda1 = diag(lamda1);
        [lamda1,ind] = sort(lamda1);
%         vec1j(jj,:,:) = vec1(:,ind);
        vNE = 0;
        for nn = 1:4
            if lamda1(nn) ~= 0
               vNE = vNE - sum(lamda1(nn)*log2(lamda1(nn)));
%              vNE = vNE - sum(lamda1(nn)*log(lamda1(nn)));
            end
        end
%         lamda1j(jj,:) = lamda1(:);
        vNEj(jj) = vNE;
        
    end
    
    Concj_ksoc(:,jk) = Concj1(:);
    vNEj_ksoc(:,jk) = vNEj(:);
    
end

% Concurrence, Fig. 5(a) and 8(a)
Concj2_ksoc = zeros(Nomega,length(ksocarray));
for jk = 1:length(ksocarray)
    for jj=1:Nomega
        Concj2_ksoc(jj,jk) = real(Concj_ksoc(jj,jk));
        if real(Concj_ksoc(jj,jk))<0
           Concj2_ksoc(jj,jk) = NaN;
        end
    end
end

Concj3_ksoc = zeros(Nomega,length(ksocarray));
for jk = 1:length(ksocarray)
    for jj=1:Nomega
        Concj3_ksoc(jj,jk) = real(Concj2_ksoc(jj,jk));
        if real(Concj2_ksoc(jj,jk))<1e-5
           Concj3_ksoc(jj,jk) = NaN;
        end
    end
end

figure
pcolor(Omegaj,ksocarray,real(Concj3_ksoc)')
set(gca,'colorscale','log')
set(gca,'FontSize',32,'FontName','Times New Roman')
shading flat
colorbar

% von Neuman entropy, Fig. 5(b) and 8(b)
figure
pcolor(Omegaj,ksocarray,real(vNEj_ksoc)')
set(gca,'FontSize',32,'FontName','Times New Roman')
shading flat
colorbar
