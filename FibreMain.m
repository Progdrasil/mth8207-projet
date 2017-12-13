% FibreMain
clear all
clc

%% Constantes
ng = 1.457420;
nc = 1.462420;
rho = 8.335e-6;
lambda = 0.6328e-6;
k=2*pi/lambda;
V=k*rho*sqrt(nc^2-ng^2);

%% Ce que tu veux modifier
%nEls, qui est le nombre d'éléments
%sinon les solutions des EF sont dans phi et betas
%Les trois premier betas sont les solutions physique. Leur fonction d'onde
%correspondantes sont phi(:,post(1:3)). Chaque colonne de phi correspond a
%une fonction et post te donne la position de la fonction. Donc si tu veux
%la fonction d'onde associé au deuxième beta, ca sera la colonne post(2).
%Sinon, c'est normal si l'amplitude entre les solutions réelles et d'ÉF
%sont différentes. Il faut normaliser. Tu divises par le max ou le min la
%solution réelle et *-1 dépendament du cas et tu obtiendras la même
%expression. c'est normal.

%Sinon, à la toute fin, tu as les valeurs de uEF, duEF et xEF correspondant
%à la solution des éléments finis.


%% Parametres
m = 2;
xMin=0;
xMax=rho*m;
nEls=m*20; %Doit être le même multiple que xMax pour qu'un node soit à rho
xn = 0;
L=0;

%% Créer le mesh
[nN,nodes,connect,nB,bEls,bPts]=meshEF(xMin,xMax,nEls);

%% degrees of freedom

pDeg=zeros(nEls,1);
pDeg(:,1)=2;     %set polynomial degree for each element
pType=zeros(nEls,1);
pType(:,1)=2;   %set element type: 1=Lagrangian, 2=hierarchical

[elDof,dFreedom]=dof(nEls,pDeg,connect);


%% Calculer les matrices du problème
%Calculer la matrice de gauche d'éléments finis
[A] = main1(nEls, xMin, xMax, xn, L);
A(end,:)=0.;A(:,end)=0.;A(end,end)=1.;
[As] = A(1:end-1,1:end-1);
%Calculer la matrice de droite d'éléments finis
[B] = main2(nEls, xMin, xMax, xn);
B(end,:)=0.;B(:,end)=0.;
[Bs] = B(1:end-1,1:end-1);

%% Résoudre le problème de valeurs propres
[phi beta2] = eig(As, Bs);
beta = sqrt(beta2);

[sy sx] = size(beta2);
for n = 1:sx
beta2d(n) = beta2(n,n);
end
[beta2s, post] = sort(beta2d,'descend');
betas = sqrt(beta2s);

%% Obtenir un plot et les valeurs de u, du et x
[uEF, duEF, xEF] = postProcFibre(nEls-1,nodes,connect,elDof,dFreedom,pDeg,pType,phi(:,post(3)));

%% Obtenir la solution exacte
[r phir]=SolutionExacte(400,xMax); %(Numbre de points pour le plot,jusqu'à où va x)

