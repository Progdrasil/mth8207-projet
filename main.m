clear all;
close all;
clc;

%% Valeurs d'imortance
%%% Variables
nEls = 8;
pDegFM = 1;
VectPropre = 3;
i=1;

%%% Constantes
ng = 1.457420;
nc = 1.462420;
rho = 8.335e-6;
lambda = 0.6328e-6;
k=2*pi/lambda;
V=k*rho*sqrt(nc^2-ng^2);

%%% Parametres
m = 2;
xMin=0;
xMax=rho*m;

%% Debut du tableaux
fprintf('h \t\t ||e||L2 \t alpha_L2 \t ||e||H1 \t alpha_H1\n')

%% Calculs
[r, phir betar, betas, post, phi, uEF, duEF, xEF]=FibreMain(xMin, xMax, nEls, pDegFM, VectPropre);


%% Analyse taux de convergence
% estimation de l'erreur en L2
normeEL2(i) = L2(uEF, xEF, phir(:,VectPropre), r);

% estimation de l'erreur en H1
normeEH1(i) = H1(uEF, xEF, 1/nEls(i), phir(:,VectPropre), r);
if i > 1
	% estimation taux de convergence en L2
	TauxEL2(i) = tauxConv(normeEL2(i), normeEL2(i-1), 1/nEls(i), 1/nEls(i-1));

	% estimation taux de covergence en H1
	TauxEH1(i) = tauxConv(normeEH1(i), normeEH1(i-1), 1/nEls(i), 1/nEls(i-1));
else
	TauxEL2(i) = 0;
	TauxEH1(i) = 0;
end

fprintf('1/%4d \t %3.2e \t %5.4f \t %3.2e \t %5.4f\n', nEls(i), normeEL2(i), TauxEL2(i), normeEH1(i), TauxEH1(i));
