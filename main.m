% clc
clear all
close all
%*******************main code for the 1D finite element solver*************
h = 2.^(2:13);
xMin=-1;
xMax=1;
p = 1;

fprintf('h \t\t ||u||L2 \t ||e||L2 \t log||e||L2 \t ||u||H1 \t ||e||H1 \t log||e||H1\n')
for i = 1:length(h)
%**************************build mesh**************************************
	nEls=(xMax - xMin) * h(i);
	[nN,nodes,connect,nB,bEls,bPts]=mesh(xMin,xMax,nEls);

%*******************degrees of freedom*************************************

	pDeg=zeros(nEls,1);
	pDeg(:,1)=p;     %set polynomial degree for each element
	%pDeg=[1,2,2,1]';   %could set each element individually
	pType=zeros(nEls,1);
	pType(:,1)=2;   %set element type: 1=Lagrangian, 2=hierarchical

	[elDof,dFreedom]=dof(nEls,pDeg,connect);

%*****************build element matrices k and vectors f*******************

	%get quadrature points and weights
	[xiQ,wQ]=gQuad;

	%build k and f for each element; sum to find K and F
	%set k,c,b,f:
	syms x
	k=@(x)-1;
	c=@(x)0;
	b=@(x)0;
	f=@(x)2;
	%f=@(x)(heaviside(x-.3)-heaviside(x-.6))*-1;
	[K,F]=element(nEls,nodes,connect,xiQ,wQ,pDeg,pType,elDof,dFreedom,k,c,b,f);

%*****************add boundary conditions and solve************************

	%set boundary condition type (set values within boundaryC.m)
	bType=zeros(2,1);   %boundary condition types; use 1 for Dirichlet,
	bType(1)=1;         %2 for Neumann, 3 for Robin
	bType(2)=1;

	[K,F]=boundaryC(nB,bEls,bPts,bType,dFreedom,K,F);
	u=K\F;

%**************post-processing*********************************************
	% calcul de lequation reel
	xt = xMin:1/(h(i) * p):xMax;
	x_true = xMin:0.001:xMax;
	xt = xt';
	ut = xt .^ 2;
	u_true = x_true .^ 2;

	% estimation de l'erreur en L2
	normeUL2(i) = L2(u, xt);
	[normeEL2(i), TauxEL2(i)] = tauxConv(ut, u, xt, 1/h(i), p, 'l2');

	% estimation de l'erreur en H1
	normeUH1(i) = H1(ut, xt, 1/h(i));
	[normeEH1(i), TauxEH1(i)] = tauxConv(ut, u, xt, 1/h(i), p, 'h1');

	fprintf('1/%4d \t %8.7f \t %3.2e \t %5.4f \t\t %8.7f \t %3.2e \t %5.4f\n', h(i), normeUL2(i), normeEL2(i), TauxEL2(i), normeUH1(i), normeEH1(i), TauxEH1(i));
end

figure
title('taux de convergence L_2')
hold on
loglog(log(h), TauxEL2)
loglog(log(h), TauxEH1)
hold off
legend('L2', 'H1', 'location', 'best')
