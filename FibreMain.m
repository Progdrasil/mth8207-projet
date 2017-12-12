% FibreMain
clear all
clc

%Constantes
ng = 1.457420;
nc = 1.462420;
rho = 8.335e-6;
lambda = 0.6328e-6;
k=2*pi/lambda;
V=k*rho*sqrt(nc^2-ng^2);

%Parametres
xMin=0;
xMax=rho*2;
nEls=200;

%Calculer la matrice de gauche d'éléments finis
[A] = main1(nEls, xMin, xMax);
A(end,:)=0.;A(:,end)=0.;A(end,end)=1.;
[As] = A(1:end-1,1:end-1);
%Calculer la matrice de droite d'éléments finis
[B] = main2(nEls, xMin, xMax);
B(end,:)=0.;B(:,end)=0.;
[Bs] = B(1:end-1,1:end-1);

%Matrice de notre problème de valeurs prores
%C = B\A;

%Conditons frontières?
%C(1,:) = 0;C(end,:) = 0;

[phi beta2] = eig(As, Bs);
beta = sqrt(beta2);

for n = 1:400
beta2d(n) = beta2(n,n);
end
[beta2s, post] = sort(beta2d,'descend');
betas = sqrt(beta2s);




figure
plot(linspace(0,rho,nEls*2),beta(388,388)*phi(:,388))
