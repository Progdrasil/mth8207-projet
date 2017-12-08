function [K,F]=element(nEls,nodes,connect,xiQ,wQ,pDeg,pType,elDof,dFreedom,k,c,b,f)
%build element matrices k & f, then K & F

m=max(max(dFreedom));

K=zeros(m);
F=zeros(m,1);
for ne=1:nEls
   x1=nodes(connect(ne,1));
   x2=nodes(connect(ne,2));
   h=x2-x1;
   jac=h/2;
   nQ=3;  %order of gaussian quadrature
   kEl=zeros(elDof(ne));
   fEl=zeros(elDof(ne),1);
   for l=1:nQ
       xi=xiQ(l,nQ);
       wq=wQ(l,nQ);
       xq=h/2*xi+(x1+x2)/2;   %global coordinate of xi
       [N,dN]=shape(xi,ne,pDeg,pType);
       dN=dN/jac;
       [kQ,cQ,bQ,fQ]=evals(xq,k,c,b,f);
       for i=1:elDof(ne)
           fEl(i)=fEl(i)+wq*fQ*N(i)*jac;   %sum element f  
           for j=1:elDof(ne)
               kEl(i,j)=kEl(i,j)+wq*kQ*dN(j)*dN(i)*jac;
               kEl(i,j)=kEl(i,j)+wq*cQ*dN(j)*N(i)*jac;
               kEl(i,j)=kEl(i,j)+wq*bQ*N(j)*N(i)*jac;
           end
       end
   end

   for i=1:elDof(ne)    %build global matrices K and F
       ci=dFreedom(ne,i);
       F(ci)=F(ci)+fEl(i);
       for j=1:elDof(ne)
           cj=dFreedom(ne,j);
           K(ci,cj)=K(ci,cj)+kEl(i,j);
       end
   end

end
   
   
   
   
   