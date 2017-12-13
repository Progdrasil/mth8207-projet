function [r phir]=SolutionExacte(nEls,mfibre)
%Solution exacte

%Constantes
nEls=50*4;
ng = 1.457420;
nc = ng+0.005;
rho = 8.335e-6;
lambda = 0.6328e-6;
k=2*pi/lambda;
V=k*rho*sqrt(nc^2-ng^2);




U(1) = 2.1845; U(2) = 4.9966; U(3) = 7.7642;
W = sqrt(V.^2-U.^2);
betar = sqrt(k^2*nc^2-(U./rho).^2);


figure
for n=1:3
    phi1=@(x)besselj(0,U(n).*x./rho)./besselj(0,U(n));
    phi2=@(x)besselk(0,W(n).*x./rho)./besselk(0,W(n));
    
    r = 0:rho/(nEls-1):mfibre;
    nb = 1;
    for x = 0:rho/(nEls-1):mfibre
        if x < rho
            phir(nb,n) = phi1(x);
        else
            phir(nb,n) = phi2(x);
        end
        nb = nb+1;
    end
%     Nlm = -rho/4*betar(n)/k*sqrt(8.854187817e-12/1.2566370614e-6);
%     Nlm = Nlm*2*pi*besselj(-1,U(n))*besselj(1,U(n))/(besselj(0,U(n)).^2);
%     NLM = Nlm*V.^2./
    phirN(:,n) = phir(:,n)./max(abs(phir(:,n)));
    if phirN(1,n) <= 0
        phirN(:,n) = -phirN(:,n)
    end
    hold on
    plot(r,phir(:,n))
end
hold off
betar

end