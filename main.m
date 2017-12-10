clc
clear all
close all
%*******************main code for the 1D finite element solver*************
h = 2.^(2:13);
xMin=0;
xMax=1;
p = 1;
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
	kt = eps;
	ut = xt .^ 2;
	u_true = x_true .^ 2;

	%plot the solution u and its derivative
	% figure(i)
	% postProc(nEls,nodes,connect,elDof,dFreedom,pDeg,pType,u);
	% hold on
	% plot(xt, u)
	% plot(x_true, u_true)
	% hold off

	% estimation de l'erreur
	h_max = 1/h(i);
	C = L2(ut - u, xt);

	normeE(i) = C * h_max^p;
	TauxE(i) = -1 * log(C) - p*log(1/h_max);
	fprintf('1/%d, \t %5.4e, \t %5.4f\n', h(i), normeE(i), TauxE(i));
end

figure
title('taux de convergence')
loglog(log(h), TauxE)
