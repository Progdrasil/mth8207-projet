function [K,F]=boundaryC(nB,bEls,bPts,bType,dFreedom,K,F)

m=max(max(dFreedom));
Kbc=K;
Fbc=F;
%set dirichlet boundary values
u0(:,1)=zeros(m,1);  
u0(:,2)=zeros(m,1); %set different boundary points to different 
u0(1,1) = 0;
u0(m,2) = 0;

%set Neumann boundary values
gamma=zeros(nB,1);
alpha=ones(nB,1);
k=1;

%set Robin boundary values      ROBIN NOT FULLY IMPLEMENTED YET
%betaR=1;
%gammaR=1;

for nb=1:nB
    be=bEls(nb);
    bp=bPts(nb);
    bt=bType(nb);
    c=dFreedom(be,bp);
    if bt==1            %Dirichlet
        Fbc=Fbc-K*u0(:,nb);
        Fbc(c)=u0(c,nb);
        Kbc(c,:)=0;
        Kbc(:,c)=0;
        Kbc(c,c)=1;
    end
    
    if bt==2            %Neumann
        Fbc(c)=Fbc(c)-k*gamma(nb)/alpha(nb);
    end
   
   % if bt==3            %Robin
   %     Kbc(c,c)=Kbc(c,c)+betaR;
   %     Fbc(c)=Fbc(c)+gammaR;
   % end
     
end
K=Kbc;
F=Fbc;
