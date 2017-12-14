% clc
clear all
close all
%*******************main code for the 1D finite element solver*************
h(1) = 4;
xMin=0;
xMax=1;
p = 1;

fprintf('h \t\t ||e||L2 \t alpha_L2 \t ||e||H1 \t alpha_H1\n')
for i = 1:10
	if i > 1
		h(i) = h(i-1) * 2;
	end
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
	figure
	hold on
	plot(x_true, u_true)
	plot(xt, u)
	hold off

	% estimation de l'erreur en L2
	normeEL2(i) = L2(u - ut, xt);

	% estimation de l'erreur en H1
	normeEH1(i) = H1(u - ut, xt, 1/h(i));
	if i > 1
		% estimation taux de convergence en L2
		TauxEL2(i) = tauxConv(normeEL2(i), normeEL2(i-1), 1/h(i), 1/h(i-1));

		% estimation taux de covergence en H1
		TauxEH1(i) = tauxConv(normeEH1(i), normeEH1(i-1), 1/h(i), 1/h(i-1));
	else
		TauxEL2(i) = 0;
		TauxEH1(i) = 0;
	end

	fprintf('1/%4d \t %3.2e \t %5.4f \t %3.2e \t %5.4f\n', h(i), normeEL2(i), TauxEL2(i), normeEH1(i), TauxEH1(i));
end

figure
title('taux de convergence L_2')
hold on
plot(h(2:end), TauxEL2(2:end))
plot(h(2:end), TauxEH1(2:end))
hold off
legend('L2', 'H1', 'location', 'best')
