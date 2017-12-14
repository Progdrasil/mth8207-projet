function [K]=main1(nEls, xMin, xMax, xn, L, p)
clc
%*******************main code for the 1D finite element solver*************
%		p = vecteur de taille 2
%			P(1) = degree polynomial
%			P(2) = type d'element (1: lagrangian, 2: hierarchical)
%**************************build mesh**************************************

% alpha=1;
% kkk=10*epsilon/alpha;
% uee=@(x)(exp((x-1)./kkk)-1)./(exp(-1/kkk)-1);
% xx=0:1/nEls:1;xx=xx';

% xMin=0;
% xMax=1.2e-6;
%nEls=100;
[nN,nodes,connect,nB,bEls,bPts]=meshEF(xMin,xMax,nEls);
%nodes=uee(1-xx);

%*******************degrees of freedom*************************************

pDeg=zeros(nEls,1);
pDeg(:,1)=p(1);     %set polynomial degree for each element
%pDeg=[1,2,2,1]';   %could set each element individually
pType=zeros(nEls,1);
pType(:,1)=p(2);   %set element type: 1=Lagrangian, 2=hierarchical

[elDof,dFreedom]=dof(nEls,pDeg,connect);

%*****************build element matrices k and vectors f*******************

%get quadrature points and weights
[xiQ,wQ]=gQuad;

%build k and f for each element; sum to find K and F
%set k,c,b,f:
syms x
k=@(x)-x.^(xn);
c=@(x)1*x.^(xn-1);
b=@(x)(((1.462420-0.005.*heaviside(x-8.335e-6)).*2.*pi./0.6328e-6).^2).*x.^(xn)+L.^2.*x.^(xn-2);
f=@(x)0;
%f=@(x)(heaviside(x-.5)-heaviside(x-.6))*-1;
[K,F]=element(nEls,nodes,connect,xiQ,wQ,pDeg,pType,elDof,dFreedom,k,c,b,f);

% %*****************add boundary conditions and solve************************
%
% %set boundary condition type (set values within boundaryC.m)
% bType=zeros(2,1);   %boundary condition types; use 1 for Dirichlet,
% bType(1)=1;         %2 for Neumann, 3 for Robin
% bType(2)=2;
%
% [K2,F2]=boundaryC(nB,bEls,bPts,bType,dFreedom,K,F);
% u=K2\F2;
%
% %**************post-processing*********************************************
% %plot the solution u and its derivative
% figure(1)
% [du]=postProc(nEls,nodes,connect,elDof,dFreedom,pDeg,pType,u);
% for i=1:nEls
% uu(i)=(u(i)+u(i+1))/2;
% uu=uu';
% dFreedom=dFreedom(end,1);
% end
% end
