clear all;
close all;
clc;

%% Valeurs d'imortance
%%% Variables
nEls = 2.^(2:8);
p = 1:3;
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
h = (xMax - xMin) ./ nEls;

%%% Solution exacte
[r phir betar dphir]=SolutionExacte(400,xMax); %(Numbre de points pour le plot,jusqu'� o� va x)

for j = 1:length(p)
	%% Reinit les vecteurs
	normeEL2 = [];
	normeEH1 = [];
	EBeta 	= [];
	TauxEL2 = [];
	TauxEH1 = [];
	TauxBeta = [];
	%% Debut des tableaux
	fprintf('\nSolutions pour P = %1d\n', p(j))
	fprintf('h \t\t ||e||L2 \t alpha_L2 \t ||e||H1 \t alpha_H1 \t ||e||beta \t alpha_beta\n')

	for i = 1:length(nEls)
		%% Calculs
		if (nEls(i) == 64)
			fig = true;
		else
			fig = false;
		end
		%%% Element finis
		[betas, post, phi, uEF, duEF, xEF]=FibreMain(xMin, xMax, nEls(i), p(j), VectPropre, phir(:, VectPropre), dphir(:, VectPropre), r, fig);

		%% Analyse taux de convergence
		% estimation de l'erreur en L2
		normeEL2(i) = L2(uEF, xEF, phir(:,VectPropre), r, nEls(i), xMax, xMin);

		% estimation de l'erreur en H1
		normeEH1(i) = H1(uEF, duEF, xEF, phir(:,VectPropre), dphir(:,VectPropre), r, nEls(i), xMax, xMin);

		% erreur quantiter d'interets
		EBeta(i) = abs(betar(VectPropre) - betas(VectPropre));
		if i > 1
			% estimation taux de convergence en L2
			TauxEL2(i) = tauxConv(normeEL2(i), normeEL2(i-1), h(i), h(i-1));

			% estimation taux de covergence en H1
			TauxEH1(i) = tauxConv(normeEH1(i), normeEH1(i-1), h(i), h(i-1));

			% estimation taux de covergence en H1
			TauxBeta(i) = tauxConv(EBeta(i), EBeta(i-1), 1/nEls(i), h(i-1));
		else
			TauxEL2(i) = 0;
			TauxEH1(i) = 0;
			TauxBeta(i) = 0;
		end

		fprintf('1/%4d \t %3.2e \t %5.4f \t %3.2e \t %5.4f \t %3.2e \t %5.4f\n', nEls(i), normeEL2(i), TauxEL2(i), normeEH1(i), TauxEH1(i), EBeta(i), TauxBeta(i));
	end
	figErr = figure;
	subplot(2, 1, 1)
	title(strcat('Analyse erreurs avec p = ', num2str(p(j))))
	hold on
	plot(1./h, TauxEL2)
	plot(1./h, TauxEH1)
	hold off
	legend('L2', 'H1', 'location', 'best')
	ylabel('\alpha(h)')
	xlabel('1/h')

	subplot(2, 1, 2)
	hold on
	loglog(1./h, normeEL2)
	loglog(1./h, normeEH1)
	hold off
	legend('L2', 'H1', 'location', 'best')
	ylabel('log ||e||')
	xlabel('log 1/h')

	figBeta = figure;
	subplot(2, 1, 1)
	title(strcat('Analyse Vecteurs Propres avec p = ', num2str(p(j))))
	plot(1./h, TauxBeta)
	ylabel('\alpha(h)')
	xlabel('1/h')

	subplot(2, 1, 2)
	loglog(1./h, EBeta)
	ylabel('log ||e_{\beta}||')
	xlabel('log 1/h')

	saveImg(figErr, strcat('err_p', num2str(p(j))));
	saveImg(figBeta, strcat('beta_p', num2str(p(j))));
end
