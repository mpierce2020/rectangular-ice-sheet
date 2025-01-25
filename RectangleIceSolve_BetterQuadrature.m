% RectangleIceSolve
% clear; close all
%% Read in plate parameters from plateeig.f
FilePath = "C:\Users\MaxPierce\Documents\Code\3dw\results\clean\corr_taper\reftran\FromSupercloud\rectangleice\io\";
Input = readtable(FilePath + "default.inp",'FileType','text');
nTrunc = Input{1,1};
nu = Input{3,1};
aob = Input{5,1};
pTrunc = Input{7,1};
checks = 0;
modebymode = 0;

% Use nondimensionalization in SquareIce sheet of YablonovitchSheet Excel
K = 1;
a = 18.190673904035528;
b = a/aob;
theta0 = pi/2;
delta = 0.052432165;
% D = 0.745604561;
% h = 0.203221994;
Gamma = 0.5;
AbsErr = 1e-2;
endrho = 50;

% Calculating k_ice and D for a given Gamma
D_find = @(D) D*kice_find(K,D,delta).^4 - Gamma;
D = fzero(D_find,1);
kice = kice_find(K,D,delta);
Gamma = D*kice^4;

%% Read in matrices of alpha from plateeig.f
A = readmatrix(FilePath + "fort.69",'FileType','text');

alphaSS = zeros(nTrunc + 1, nTrunc + 1, pTrunc);
alphaSA = zeros(nTrunc + 1, nTrunc + 1, pTrunc);
alphaAS = zeros(nTrunc + 1, nTrunc + 1, pTrunc);
alphaAA = zeros(nTrunc + 1, nTrunc + 1, pTrunc);

for p = 1:pTrunc
    alphaSS(:,:,p) = reshape(A((1:(nTrunc+1)^2) + (p-1)*(nTrunc+1)^2,4), nTrunc+1, nTrunc+1);
    % alphaSS(:,:,p) = alphaSS(:,:,p)';
    alphaSA(:,:,p) = reshape(A((1:(nTrunc+1)^2) + (p-1)*(nTrunc+1)^2,5), nTrunc+1, nTrunc+1);
    % alphaSA(:,:,p) = alphaSA(:,:,p)';
    alphaAS(:,:,p) = reshape(A((1:(nTrunc+1)^2) + (p-1)*(nTrunc+1)^2,6), nTrunc+1, nTrunc+1);
    % alphaAS(:,:,p) = alphaAS(:,:,p)';
    alphaAA(:,:,p) = reshape(A((1:(nTrunc+1)^2) + (p-1)*(nTrunc+1)^2,7), nTrunc+1, nTrunc+1);
    % alphaAA(:,:,p) = alphaAA(:,:,p)';
end

Lnondim = readmatrix(FilePath + "default.res",'FileType','text');
Ldim = Lnondim/a^4;
lambda_SS = Ldim(:,2);
lambda_SA = Ldim(:,3);
lambda_AS = Ldim(:,4);
lambda_AA = Ldim(:,5);

%% Read in in vacuo wavenumbers
kRead = readmatrix(FilePath + "fort.70",'FileType','text');
k = kRead(:,2);

%% Checks
if checks == 1
    dx = a/50;
    x = -a:dx:a;
    y = -b:dx:b;
    [X,Y] = meshgrid(x,y);
    % Checking that in vacuo modes are orthogonal
    ESS = zeros(pTrunc,pTrunc);
    ESA = zeros(pTrunc,pTrunc);
    EAS = zeros(pTrunc,pTrunc);
    EAA = zeros(pTrunc,pTrunc);
    for ii = 1:pTrunc
        for jj = 1:pTrunc
            ESS(ii,jj) = 1/(a*b)*trapz(y,trapz(x,CapW(nTrunc, k, alphaSS, a, b, X, Y, ii, 0, 0).*CapW(nTrunc, k, alphaSS, a, b, X, Y, jj, 0, 0),1),2);
            ESA(ii,jj) = 1/(a*b)*trapz(y,trapz(x,CapW(nTrunc, k, alphaSA, a, b, X, Y, ii, 0, 1).*CapW(nTrunc, k, alphaSA, a, b, X, Y, jj, 0, 1),1),2);
            EAS(ii,jj) = 1/(a*b)*trapz(y,trapz(x,CapW(nTrunc, k, alphaAS, a, b, X, Y, ii, 1, 0).*CapW(nTrunc, k, alphaAS, a, b, X, Y, jj, 1, 0),1),2);
            EAA(ii,jj) = 1/(a*b)*trapz(y,trapz(x,CapW(nTrunc, k, alphaAA, a, b, X, Y, ii, 1, 1).*CapW(nTrunc, k, alphaAA, a, b, X, Y, jj, 1, 1),1),2);
        end
    end
   
    K = 1;
    a = 1;
    b = 1;
    dchi = 0.01;
    chi = (0 + dchi):dchi:pi/2;
    mode1 = 1;
    mode2 = 2;
    Int1 = trapz(chi,wScriptEven(K*a*cos(chi),mode1,k) .* wScriptEven(K*b*sin(chi),mode1,k) .* wScriptEven(-K*a*cos(chi),mode2,k) .* wScriptEven(-K*b*sin(chi),mode2,k));
    Int1func = @() Int1;
    t1 = timeit(Int1func);
    Int1_integrand = @(chi) wScriptEven(K*a*cos(chi),mode1,k) .* wScriptEven(K*b*sin(chi),mode1,k) .* wScriptEven(-K*a*cos(chi),mode2,k) .* wScriptEven(-K*b*sin(chi),mode2,k);
    Int1gauss = integral(Int1_integrand,0,pi/2);
    Int1gaussfunc = @() Int1gauss;
    t2 = timeit(Int1gaussfunc);
    % Int1 = trapz(chi,wScriptOdd(K*a*cos(chi),mode1,k) .* wScriptOdd(K*b*sin(chi),mode1,k) .* wScriptOdd(-K*a*cos(chi),mode2,k) .* wScriptOdd(-K*b*sin(chi),mode2,k));
    chi = 0.2;
    drho = 0.005;
    rho1 = drho:drho:1-drho;
    rho2 = (2 + drho):drho:50;
    Int2 = trapz(rho1,rho1.^2./(rho1-1).*wScriptEven(rho1*cos(chi),mode1,k) .* wScriptEven(rho1*sin(chi),mode1,k) .* wScriptEven(-rho1*cos(chi),mode2,k) .* wScriptEven(-rho1*sin(chi),mode2,k) + ...
        (2-rho1).^2./((2-rho1)-1).*wScriptEven((2-rho1)*cos(chi),mode1,k) .* wScriptEven((2-rho1)*sin(chi),mode1,k) .* wScriptEven(-(2-rho1)*cos(chi),mode2,k) .* wScriptEven(-(2-rho1)*sin(chi),mode2,k)) + ...
        trapz(rho2,rho2.^2./(rho2-1).*wScriptEven(rho2*cos(chi),mode1,k) .* wScriptEven(rho2*sin(chi),mode1,k) .* wScriptEven(-rho2*cos(chi),mode2,k) .* wScriptEven(-rho2*sin(chi),mode2,k));
    Int_part2 = trapz(rho2,rho2.^2./(rho2-1).*wScriptEven(rho2*cos(chi),mode1,k) .* wScriptEven(rho2*sin(chi),mode1,k) .* wScriptEven(-rho2*cos(chi),mode2,k) .* wScriptEven(-rho2*sin(chi),mode2,k));
    figure
    plot(rho2,rho2.^2./(rho2-1).*wScriptEven(rho2*cos(chi),mode1,k) .* wScriptEven(rho2*sin(chi),mode1,k) .* wScriptEven(-rho2*cos(chi),mode2,k) .* wScriptEven(-rho2*sin(chi),mode2,k))
    % tst = wScriptEven(2*cos(chi),mode1,k) .* wScriptEven(2*sin(chi),mode1,k) .* wScriptEven(-2*cos(chi),mode2,k) .* wScriptEven(-2*sin(chi),mode2,k)
    % tst2 = rho2(1).^2/(rho2(1)-1)
    return
end

%% Forming Q_ij

% Q_SS  = zeros(pTrunc,pTrunc);
% Q_SA  = zeros(pTrunc,pTrunc);
% Q_AS  = zeros(pTrunc,pTrunc);
% Q_AA  = zeros(pTrunc,pTrunc);

% only looping through upper diagonal since Q_ij is known to be symmetric
for ii = 36 %1:pTrunc
    for jj = 36 %ii:pTrunc
        % SS
        QPV_SSintegrand1 = @(chigrid1,rhogrid1) rhogrid1.^2./(rhogrid1 - 1).*CapWScript(nTrunc, k, alphaSS, a, b, K*rhogrid1.*cos(chigrid1), K*rhogrid1.*sin(chigrid1), ii, 0, 0).*CapWScript(nTrunc, k, alphaSS, a, b, -K*rhogrid1.*cos(chigrid1), -K*rhogrid1.*sin(chigrid1), jj, 0, 0) + ...
            (2 - rhogrid1).^2./((2 - rhogrid1) - 1).*CapWScript(nTrunc, k, alphaSS, a, b, K*(2 - rhogrid1).*cos(chigrid1), K*(2 - rhogrid1).*sin(chigrid1), ii, 0, 0).*CapWScript(nTrunc, k, alphaSS, a, b, -K*(2 - rhogrid1).*cos(chigrid1), -K*(2 - rhogrid1).*sin(chigrid1), jj, 0, 0);
        QPV_SSintegrand2 = @(chigrid2, rhogrid2) rhogrid2.^2./(rhogrid2 - 1).*CapWScript(nTrunc, k, alphaSS, a, b, K*rhogrid2.*cos(chigrid2), K*rhogrid2.*sin(chigrid2), ii, 0, 0).*CapWScript(nTrunc, k, alphaSS, a, b, -K*rhogrid2.*cos(chigrid2), -K*rhogrid2.*sin(chigrid2), jj, 0, 0);
        QPV_SS_gauss = integral2(QPV_SSintegrand1, 0, pi/2, 0, 1,'Method','iterated','AbsTol',AbsErr) + integral2(QPV_SSintegrand2, 0, pi/2, 2, endrho,'Method','iterated','AbsTol',AbsErr);

        QOther_SSintegrand = @(chi) CapWScript(nTrunc, k, alphaSS, a, b, K*cos(chi),K*sin(chi), ii, 0, 0).*CapWScript(nTrunc, k, alphaSS, a, b, -K*cos(chi),-K*sin(chi), jj, 0, 0);
        QOther_SS_gauss = integral(QOther_SSintegrand,0,pi/2,'AbsTol',AbsErr);

        Q_SS(ii,jj) = K^2*a*b/pi*(1/pi*QPV_SS_gauss + 1i*QOther_SS_gauss);

        % SA
        QPV_SAintegrand1 = @(chigrid1,rhogrid1) rhogrid1.^2./(rhogrid1 - 1).*CapWScript(nTrunc, k, alphaSA, a, b, K*rhogrid1.*cos(chigrid1), K*rhogrid1.*sin(chigrid1), ii, 0, 1).*CapWScript(nTrunc, k, alphaSA, a, b, -K*rhogrid1.*cos(chigrid1), -K*rhogrid1.*sin(chigrid1), jj, 0, 1) + ...
            (2 - rhogrid1).^2./((2 - rhogrid1) - 1).*CapWScript(nTrunc, k, alphaSA, a, b, K*(2 - rhogrid1).*cos(chigrid1), K*(2 - rhogrid1).*sin(chigrid1), ii, 0, 1).*CapWScript(nTrunc, k, alphaSA, a, b, -K*(2 - rhogrid1).*cos(chigrid1), -K*(2 - rhogrid1).*sin(chigrid1), jj, 0, 1);
        QPV_SAintegrand2 = @(chigrid2,rhogrid2) rhogrid2.^2./(rhogrid2 - 1).*CapWScript(nTrunc, k, alphaSA, a, b, K*rhogrid2.*cos(chigrid2), K*rhogrid2.*sin(chigrid2), ii, 0, 1).*CapWScript(nTrunc, k, alphaSA, a, b, -K*rhogrid2.*cos(chigrid2), -K*rhogrid2.*sin(chigrid2), jj, 0, 1);
        QPV_SA_gauss = integral2(QPV_SAintegrand1, 0, pi/2, 0, 1,'Method','iterated','AbsTol',AbsErr) + integral2(QPV_SAintegrand2, 0, pi/2, 2, endrho,'Method','iterated','AbsTol',AbsErr);

        QOther_SAintegrand = @(chi) CapWScript(nTrunc, k, alphaSA, a, b, K*cos(chi),K*sin(chi), ii, 0, 1).*CapWScript(nTrunc, k, alphaSA, a, b, -K*cos(chi),-K*sin(chi), jj, 0, 1);
        QOther_SA_gauss = integral(QOther_SAintegrand, 0, pi/2, 'AbsTol',AbsErr);

        Q_SA(ii,jj) = K^2*a*b/pi*(1/pi*QPV_SA_gauss + 1i*QOther_SA_gauss);

        % AS
        QPV_ASintegrand1 = @(chigrid1,rhogrid1) rhogrid1.^2./(rhogrid1 - 1).*CapWScript(nTrunc, k, alphaAS, a, b, K*rhogrid1.*cos(chigrid1), K*rhogrid1.*sin(chigrid1), ii, 1, 0).*CapWScript(nTrunc, k, alphaAS, a, b, -K*rhogrid1.*cos(chigrid1), -K*rhogrid1.*sin(chigrid1), jj, 1, 0) + ...
            (2 - rhogrid1).^2./((2 - rhogrid1) - 1).*CapWScript(nTrunc, k, alphaAS, a, b, K*(2 - rhogrid1).*cos(chigrid1), K*(2 - rhogrid1).*sin(chigrid1), ii, 1, 0).*CapWScript(nTrunc, k, alphaAS, a, b, -K*(2 - rhogrid1).*cos(chigrid1), -K*(2 - rhogrid1).*sin(chigrid1), jj, 1, 0);
        QPV_ASintegrand2 = @(chigrid2,rhogrid2) rhogrid2.^2./(rhogrid2 - 1).*CapWScript(nTrunc, k, alphaAS, a, b, K*rhogrid2.*cos(chigrid2), K*rhogrid2.*sin(chigrid2), ii, 1, 0).*CapWScript(nTrunc, k, alphaAS, a, b, -K*rhogrid2.*cos(chigrid2), -K*rhogrid2.*sin(chigrid2), jj, 1, 0);
        QPV_AS_gauss = integral2(QPV_ASintegrand1, 0, pi/2, 0, 1,'Method','iterated', 'AbsTol', AbsErr) + integral2(QPV_ASintegrand2, 0, pi/2, 2, endrho,'Method','iterated', 'AbsTol', AbsErr);

        QOther_ASintegrand = @(chi) CapWScript(nTrunc, k, alphaAS, a, b, K*cos(chi),K*sin(chi), ii, 1, 0).*CapWScript(nTrunc, k, alphaAS, a, b, -K*cos(chi),-K*sin(chi), jj, 1, 0);
        QOther_AS_gauss = integral(QOther_ASintegrand, 0, pi/2, 'AbsTol', AbsErr);

        Q_AS(ii,jj) = K^2*a*b/pi*(1/pi*QPV_AS_gauss + 1i*QOther_AS_gauss);

         % AA
        QPV_AAintegrand1 = @(chigrid1,rhogrid1) rhogrid1.^2./(rhogrid1 - 1).*CapWScript(nTrunc, k, alphaAA, a, b, K*rhogrid1.*cos(chigrid1), K*rhogrid1.*sin(chigrid1), ii, 1, 1).*CapWScript(nTrunc, k, alphaAA, a, b, -K*rhogrid1.*cos(chigrid1), -K*rhogrid1.*sin(chigrid1), jj, 1, 1) + ...
            (2 - rhogrid1).^2./((2 - rhogrid1) - 1).*CapWScript(nTrunc, k, alphaAA, a, b, K*(2 - rhogrid1).*cos(chigrid1), K*(2 - rhogrid1).*sin(chigrid1), ii, 1, 1).*CapWScript(nTrunc, k, alphaAA, a, b, -K*(2 - rhogrid1).*cos(chigrid1), -K*(2 - rhogrid1).*sin(chigrid1), jj, 1, 1);
        QPV_AAintegrand2 = @(chigrid2,rhogrid2) rhogrid2.^2./(rhogrid2 - 1).*CapWScript(nTrunc, k, alphaAA, a, b, K*rhogrid2.*cos(chigrid2), K*rhogrid2.*sin(chigrid2), ii, 1, 1).*CapWScript(nTrunc, k, alphaAA, a, b, -K*rhogrid2.*cos(chigrid2), -K*rhogrid2.*sin(chigrid2), jj, 1, 1);
        QPV_AA_gauss = integral2(QPV_AAintegrand1, 0, pi/2, 0, 1,'Method','iterated', 'AbsTol', AbsErr) + integral2(QPV_AAintegrand2, 0, pi/2, 2, endrho,'Method','iterated', 'AbsTol', AbsErr);

        QOther_AAintegrand = @(chi) CapWScript(nTrunc, k, alphaAA, a, b, K*cos(chi),K*sin(chi), ii, 1, 1).*CapWScript(nTrunc, k, alphaAA, a, b, -K*cos(chi),-K*sin(chi), jj, 1, 1);
        QOther_AA_gauss = integral(QOther_AAintegrand, 0, pi/2, 'AbsTol', AbsErr);

        Q_AA(ii,jj) = K^2*a*b/pi*(1/pi*QPV_AA_gauss + 1i*QOther_AA_gauss);
    end
end
% reflecting matrices over diagonal to make symmetric
% Q_SS = Q_SS + triu(Q_SS, 1).';
% Q_SA = Q_SA + triu(Q_SA, 1).';
% Q_AS = Q_AS + triu(Q_AS, 1).';
% Q_AA = Q_AA + triu(Q_AA, 1).';

R_SS = zeros(pTrunc,1);
R_SA = zeros(pTrunc,1);
R_AS = zeros(pTrunc,1);
R_AA = zeros(pTrunc,1);
if modebymode == 0
for ii = 1:pTrunc
    R_SS(ii) = CapWScript(nTrunc, k, alphaSS, a, b, -K*cos(theta0), -K*sin(theta0), ii, 0, 0);
    R_SA(ii) = CapWScript(nTrunc, k, alphaSA, a, b, -K*cos(theta0), -K*sin(theta0), ii, 0, 1);
    R_AS(ii) = CapWScript(nTrunc, k, alphaAS, a, b, -K*cos(theta0), -K*sin(theta0), ii, 1, 0);
    R_AA(ii) = CapWScript(nTrunc, k, alphaAA, a, b, -K*cos(theta0), -K*sin(theta0), ii, 1, 1);
end
elseif modebymode == 1
    R_SA(3) = 1/(K*delta - D*lambda_SA(3));
end

E_SS = zeros(pTrunc,1);
E_SA = zeros(pTrunc,1);
E_AS = zeros(pTrunc,1);
E_AA = zeros(pTrunc,1);
for ii = 1:pTrunc
    E_SS(ii) = CapE(nTrunc, alphaSS, ii);
    E_SA(ii) = CapE(nTrunc, alphaSA, ii);
    E_AS(ii) = CapE(nTrunc, alphaAS, ii);
    E_AA(ii) = CapE(nTrunc, alphaAA, ii);
end

LHS_firstterm_SS = diag(E_SS./(K*delta - D*lambda_SS(1:pTrunc)));
LHS_firstterm_SA = diag(E_SA./(K*delta - D*lambda_SA(1:pTrunc)));
LHS_firstterm_AS = diag(E_AS./(K*delta - D*lambda_AS(1:pTrunc)));
LHS_firstterm_AA = diag(E_AA./(K*delta - D*lambda_AA(1:pTrunc)));

LHS_SS = LHS_firstterm_SS - Q_SS.';
LHS_SA = LHS_firstterm_SA - Q_SA.';
LHS_AS = LHS_firstterm_AS - Q_AS.';
LHS_AA = LHS_firstterm_AA - Q_AA.';

a_SS = LHS_SS\R_SS;
a_SA = LHS_SA\R_SA;
a_AS = LHS_AS\R_AS;
a_AA = LHS_AA\R_AA;

dx = a/50;
x = -a:dx:a;
y = -b:dx:b;
xtot = -3*a:dx:3*a;
ytot = -3*b:dx:3*b;
[X,Y] = meshgrid(x,y);
[Xtot,Ytot] = meshgrid(xtot,ytot);
etaincident = exp(1i*K*(Xtot*cos(theta0) + Ytot*sin(theta0)));
% etaincident((round(size(Xtot,1)/4):round(3*size(Xtot,1)/4)), (round(size(Xtot,1)/4):round(3*size(Xtot,1)/4))) = zeros(length((round(size(Xtot,1)/4):round(3*size(Xtot,1)/4))), length((round(size(Xtot,1)/4):round(3*size(Xtot,1)/4))));
Inner = padarray(ones(length(x),length(y)),[length(x)-1,length(y)-1],0);
Outer = padarray(zeros(length(x),length(y)),[length(x)-1,length(y)-1],1);

etaplateSS = zeros(length(xtot),length(ytot));
etaplateSA = zeros(length(xtot),length(ytot));
etaplateAS = zeros(length(xtot),length(ytot));
etaplateAA = zeros(length(xtot),length(ytot));
etaplateonly = zeros(length(x),length(y));

for ii = 1:pTrunc
    etaplateSS = etaplateSS + a_SS(ii)/(K*delta - D*lambda_SS(ii))*CapW(nTrunc,k,alphaSS,a,b,Xtot,Ytot,ii,0,0);
    etaplateSA = etaplateSA + a_SA(ii)/(K*delta - D*lambda_SA(ii))*CapW(nTrunc,k,alphaSA,a,b,Xtot,Ytot,ii,0,1);
    etaplateAS = etaplateAS + a_AS(ii)/(K*delta - D*lambda_AS(ii))*CapW(nTrunc,k,alphaAS,a,b,Xtot,Ytot,ii,1,0);
    etaplateAA = etaplateAA + a_AA(ii)/(K*delta - D*lambda_AA(ii))*CapW(nTrunc,k,alphaAA,a,b,Xtot,Ytot,ii,1,1);

    etaplateonly = etaplateonly + a_SS(ii)/(K*delta - D*lambda_SS(ii))*CapW(nTrunc,k,alphaSS,a,b,X,Y,ii,0,0) + ...
        a_SA(ii)/(K*delta - D*lambda_SA(ii))*CapW(nTrunc,k,alphaSA,a,b,X,Y,ii,0,1) + ...
        a_AS(ii)/(K*delta - D*lambda_AS(ii))*CapW(nTrunc,k,alphaAS,a,b,X,Y,ii,1,0) + ...
        a_AA(ii)/(K*delta - D*lambda_AA(ii))*CapW(nTrunc,k,alphaAA,a,b,X,Y,ii,1,1);
end

etaplatesquared_spatialavg = mean(etaplateonly.*conj(etaplateonly),'all');
etaplate_spatialavg = mean(etaplateonly,'all');
b0squared_analytical_yablonovitch = 1/((1+Gamma)^2*(1+5*Gamma));

etaplate = etaplateSS + etaplateSA + etaplateAS + etaplateAA;
eta = etaplate.*Inner + etaincident.*Outer;

phs = 0;

figure(1)
surf(Xtot,Ytot,real(exp(1i*phs)*eta),'EdgeColor','none')
xlabel('x')
ylabel('y')

figure(2)
surf(Xtot,Ytot,real(exp(1i*phs)*eta),'EdgeColor','none')
xlim([-a,a])
ylim([-b,b])
xlabel('x')
ylabel('y')

% Checking using optical theorem and comparing figures with Porter paper
theta = 0:0.01:2*pi;
DiffCoeff_unscaled = zeros(1,length(theta));
DiffCoeffTheta0_unscaled = 0;
theta_fwd = mod(theta0 + pi, 2*pi);

for ii = 1:pTrunc
    DiffCoeff_unscaled = DiffCoeff_unscaled + a_SS(ii)*CapWScript(nTrunc,k,alphaSS,a,b,-K*cos(theta),-K*sin(theta),ii,0,0) + ...
        a_SA(ii)*CapWScript(nTrunc,k,alphaSA,a,b,-K*cos(theta),-K*sin(theta),ii,0,1) + ...
        a_AS(ii)*CapWScript(nTrunc,k,alphaAS,a,b,-K*cos(theta),-K*sin(theta),ii,1,0) + ...
        a_AA(ii)*CapWScript(nTrunc,k,alphaAA,a,b,-K*cos(theta),-K*sin(theta),ii,1,1);
    DiffCoeffTheta0_unscaled = DiffCoeffTheta0_unscaled + a_SS(ii)*CapWScript(nTrunc,k,alphaSS,a,b,-K*cos(theta_fwd),-K*sin(theta_fwd),ii,0,0) + ...
        a_SA(ii)*CapWScript(nTrunc,k,alphaSA,a,b,-K*cos(theta_fwd),-K*sin(theta_fwd),ii,0,1) + ...
        a_AS(ii)*CapWScript(nTrunc,k,alphaAS,a,b,-K*cos(theta_fwd),-K*sin(theta_fwd),ii,1,0) + ...
        a_AA(ii)*CapWScript(nTrunc,k,alphaAA,a,b,-K*cos(theta_fwd),-K*sin(theta_fwd),ii,1,1);
end
DiffCoeff = 1i/2*K^2*a*b*DiffCoeff_unscaled;
DiffCoeffTheta0 = 1i/2*K^2*a*b*DiffCoeffTheta0_unscaled;


figure
hold on
plot(theta,abs(DiffCoeff))
plot(theta,real(DiffCoeff))
plot(theta,imag(DiffCoeff))
xline(theta0)
xline(theta_fwd,'--')
hold off
legend('Abs', 'Real', 'Imag','$\theta_0$','$\theta_{fwd}$')
xlabel('$\theta$')

figure
plot(flip(theta-theta0),abs(DiffCoeff))
xlim([0,pi])
xlabel('$\theta$')
ylabel('$|A(\theta,0)|$')
% title("$k^Ia$=" + round(kice*a,2) + ", $\Gamma=$" + round(Gamma,2) + ", $k^Oh=$" + round(K*h,2))


ScatteringCrossSection = 1/(2*pi)*trapz(theta(2:end),abs(DiffCoeff(2:end)).^2);
OptTheorem_err = ScatteringCrossSection + real(DiffCoeffTheta0);



%% Functions
function y = wScriptEven(sigma, m, k)
% +1 to k index to account for zero indexing in Fortran and one indexing in
% Matlab
if m == 0
    y = sqrt(2)*sin(sigma)./sigma;
else 
    y = 2*sqrt(2)*sigma.^2./(sigma.^4 - k(2*m+1)^4).*(k(2*m+1)*cos(sigma)*tanh(k(2*m+1)) + sigma.*sin(sigma));
end
end

function y = wScriptOdd(sigma, m, k)
% +1 to k index to account for zero indexing in Fortran and one indexing in
% Matlab
if m == 0
    y = -1i*sqrt(6)*(sin(sigma) - sigma.*cos(sigma))./sigma.^2;
else 
    y = -1i*2*sqrt(2)*sigma.^2./(sigma.^4 - k(2*m + 2).^4).*(k(2*m + 2)*sin(sigma)*coth(k(2*m + 2)) - sigma.*cos(sigma));
end
end

function y = wEven(t, m, k)
% even eigenmodes corresponding to free-free Euler-Bernoulli beam in physical
% domain
% +1 to k index to account for zero indexing in Fortran and one indexing in
% Matlab
if m == 0
    y = sqrt(1/2);
else 
    y = 1/sqrt(2)*(cosh(k(2*m+1)*t)/cosh(k(2*m+1)) + cos(k(2*m+1)*t)/cos(k(2*m+1)));
end
end

function y = wOdd(t, m, k)
% odd eigenmodes corresponding to free-free Euler-Bernoulli beam in physical
% domain
% +1 to k index to account for zero indexing in Fortran and one indexing in
% Matlab
if m == 0
    y = t*sqrt(3/2);
else 
    y = 1/sqrt(2)*(sinh(k(2*m + 2)*t)/sinh(k(2*m + 2)) + sin(k(2*m + 2)*t)/sin(k(2*m + 2)));
end
end

% function y = CapWScript(nTrunc, k, alpha, a, b, alphadum, betadum, p, mu, nu)
% y = zeros(size(betadum,1), size(betadum,2));
% if mu == 0 && nu == 0 % SS
%     for mm = 0:(nTrunc-1)
%         for nn = 0:(nTrunc-1)
%             y = y + alpha(mm+1,nn+1,p)*wScriptEven(alphadum*a,mm,k).*wScriptEven(betadum*b,nn,k);
%         end
%     end
% elseif mu == 0 && nu == 1 % SA
%     for mm = 0:(nTrunc-1)
%         for nn = 0:(nTrunc-1)
%             y = y + alpha(mm+1,nn+1,p)*wScriptEven(alphadum*a,mm,k).*wScriptOdd(betadum*b,nn,k);
%         end
%     end
% elseif mu == 1 && nu == 0 % AS
%     for mm = 0:(nTrunc-1)
%         for nn = 0:(nTrunc-1)
%             y = y + alpha(mm+1,nn+1,p)*wScriptOdd(alphadum*a,mm,k).*wScriptEven(betadum*b,nn,k);
%         end
%     end
% elseif mu == 1 && nu == 1 % AA
%     for mm = 0:(nTrunc-1)
%         for nn = 0:(nTrunc-1)
%             y = y + alpha(mm+1,nn+1,p)*wScriptOdd(alphadum*a,mm,k).*wScriptOdd(betadum*b,nn,k);
%         end
%     end
% end
% end

function y = CapWScript(nTrunc, k, alpha, a, b, alphadum, betadum, p, mu, nu)
    % Precompute reused terms to avoid redundant calculations
    alphadum_a = alphadum * a;
    betadum_b = betadum * b;
    
    % Initialize output matrix
    y = zeros(size(betadum, 1), size(betadum, 2));
    
    % Precompute wScript results for all needed indices
    wEven_alpha = arrayfun(@(mm) wScriptEven(alphadum_a, mm, k), 0:(nTrunc-1), 'UniformOutput', false);
    wEven_beta = arrayfun(@(nn) wScriptEven(betadum_b, nn, k), 0:(nTrunc-1), 'UniformOutput', false);
    wOdd_alpha = arrayfun(@(mm) wScriptOdd(alphadum_a, mm, k), 0:(nTrunc-1), 'UniformOutput', false);
    wOdd_beta = arrayfun(@(nn) wScriptOdd(betadum_b, nn, k), 0:(nTrunc-1), 'UniformOutput', false);

    % Compute based on (mu, nu) conditions
    if mu == 0 && nu == 0 % SS
        for mm = 0:(nTrunc-1)
            for nn = 0:(nTrunc-1)
                y = y + alpha(mm+1, nn+1, p) * wEven_alpha{mm+1} .* wEven_beta{nn+1};
            end
        end
    elseif mu == 0 && nu == 1 % SA
        for mm = 0:(nTrunc-1)
            for nn = 0:(nTrunc-1)
                y = y + alpha(mm+1, nn+1, p) * wEven_alpha{mm+1} .* wOdd_beta{nn+1};
            end
        end
    elseif mu == 1 && nu == 0 % AS
        for mm = 0:(nTrunc-1)
            for nn = 0:(nTrunc-1)
                y = y + alpha(mm+1, nn+1, p) * wOdd_alpha{mm+1} .* wEven_beta{nn+1};
            end
        end
    elseif mu == 1 && nu == 1 % AA
        for mm = 0:(nTrunc-1)
            for nn = 0:(nTrunc-1)
                y = y + alpha(mm+1, nn+1, p) * wOdd_alpha{mm+1} .* wOdd_beta{nn+1};
            end
        end
    end
end

function out = CapW(nTrunc, k, alpha, a, b, x, y, p, mu, nu)
out = zeros(size(x,1), size(x,2));
if mu == 0 && nu == 0 % SS
    for mm = 0:(nTrunc-1)
        for nn = 0:(nTrunc-1)
            out = out + alpha(mm+1,nn+1,p)*wEven(x/a,mm,k).*wEven(y/b,nn,k);
        end
    end
elseif mu == 0 && nu == 1 % SA
    for mm = 0:(nTrunc-1)
        for nn = 0:(nTrunc-1)
            out = out + alpha(mm+1,nn+1,p)*wEven(x/a,mm,k).*wOdd(y/b,nn,k);
        end
    end
elseif mu == 1 && nu == 0 % AS
    for mm = 0:(nTrunc-1)
        for nn = 0:(nTrunc-1)
            out = out + alpha(mm+1,nn+1,p)*wOdd(x/a,mm,k).*wEven(y/b,nn,k);
        end
    end
elseif mu == 1 && nu == 1 % AA
    for mm = 0:(nTrunc-1)
        for nn = 0:(nTrunc-1)
            out = out + alpha(mm+1,nn+1,p)*wOdd(x/a,mm,k).*wOdd(y/b,nn,k);
        end
    end
end
end

% function out = CapW(nTrunc, k, alpha, a, b, x, y, p, mu, nu)
%     % Precompute reused terms to avoid redundant calculations
%     x_scaled = x / a;
%     y_scaled = y / b;
% 
%     % Initialize output matrix
%     out = zeros(size(x, 1), size(x, 2));
% 
%     % Precompute w functions for all required indices
%     wEven_x = cell(1, nTrunc);
%     wOdd_x = cell(1, nTrunc);
%     wEven_y = cell(1, nTrunc);
%     wOdd_y = cell(1, nTrunc);
% 
%     for mm = 0:(nTrunc-1)
%         wEven_x{mm+1} = wEven(x_scaled, mm, k);
%         wOdd_x{mm+1} = wOdd(x_scaled, mm, k);
%     end
% 
%     for nn = 0:(nTrunc-1)
%         wEven_y{nn+1} = wEven(y_scaled, nn, k);
%         wOdd_y{nn+1} = wOdd(y_scaled, nn, k);
%     end
% 
%     % Compute the output based on mu and nu conditions
%     for mm = 0:(nTrunc-1)
%         for nn = 0:(nTrunc-1)
%             alpha_value = alpha(mm+1, nn+1, p);
% 
%             if mu == 0 && nu == 0 % SS
%                 out = out + alpha_value * wEven_x{mm+1} .* wEven_y{nn+1};
%             elseif mu == 0 && nu == 1 % SA
%                 out = out + alpha_value * wEven_x{mm+1} .* wOdd_y{nn+1};
%             elseif mu == 1 && nu == 0 % AS
%                 out = out + alpha_value * wOdd_x{mm+1} .* wEven_y{nn+1};
%             elseif mu == 1 && nu == 1 % AA
%                 out = out + alpha_value * wOdd_x{mm+1} .* wOdd_y{nn+1};
%             end
%         end
%     end
% end

function y = CapE(nTrunc, alpha, p)
y = 0;
    for mm = 0:(nTrunc-1)
        for nn = 0:(nTrunc-1)
            y = y + alpha(mm+1,nn+1,p).^2;
        end
    end
end

function kice = kice_find(K,D,delta)
    kice_fnc = @(kice) (kice + D*kice.^5)./(kice*delta + 1) - K;
    kice = fzero(kice_fnc, K);
end