function main
clc
close all
%*******************main code for the 1D finite element solver*************

%**************************build mesh**************************************

xMin=0;
xMax=1;
nEls=2;
[nN,nodes,connect,nB,bEls,bPts]=mesh(xMin,xMax,nEls);

%*******************degrees of freedom*************************************

pDeg=zeros(nEls,1);
pDeg(:,1)=1;     %set polynomial degree for each element
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
k=@(x)1;
c=@(x)0;
b=@(x)0;
f=@(x)-1;
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
%plot the solution u and its derivative
figure(1)
postProc(nEls,nodes,connect,elDof,dFreedom,pDeg,pType,u);

