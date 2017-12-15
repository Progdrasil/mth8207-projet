% FibreMain
function [betas, post, phi, uEF, duEF, xEF]=FibreMain(xMin, xMax, nEls, pDegFM, VectPropre, phir, dphir, r, fig)
%% Instructions
%Inputs
% nEls : DOIT ETRE UN NOMBRE PAIRE. C'est le nombre d'�l�ments
% pDegFM : Ordre de p pour une fonction de forme lagrangienne
% VectPropre : Sort les uEF, duEF, xEF pour le i�me vecteur propre. Donc,
% si tu veux faire ton analyse pour le troisi�me vecteur et valeur propre,
% tu mets =3 pour que uEF, duEF et xEF que la fonction te donne soit els
% valeurs pour le troisi�me vecteur et valeur propre

%Outputs
% uEF, duEF, xEF : Te donne ce que tu veux. u est la fonction du ca d�riv�
% et x le vecteur x correspondant. Te le donne pour un seul vecteur propre
% correpondant � la valeur de VectPropre.
% r : Vecteur x pour la solution exacte
% phir : Solution exacte des vecteur propres. Chaque colonne est un vecteur
% diff�rent. � besoin d'�tre normalis�. Divise par son max ou min, p-e *-1
% si tu veux plot les deux solutions sur une m�me �chellle.
% betar : Solution exacte des valeurs propres. Chaque colonne est une valeur diff�rent.

% betas : Valeurs propres des �l�ments finis en ordre. C-�-d, que la
% troisi�me valeur propre correspond � betas(3).
% post : Te donne la position du vecteur propre correspondant. Si tu veux
% le vecteur propre phi associ� � betas(3), tu fais phi(:,post(3)). Tu n'en
% as pas besoin puisque tu vas utiliser les uEF, duEF, xEF


%% Constantes
% ng = 1.457420;
% nc = 1.462420;
% rho = 8.335e-6;
% lambda = 0.6328e-6;
% k=2*pi/lambda;
% V=k*rho*sqrt(nc^2-ng^2);

%% Ce que tu veux modifier
%nEls, qui est le nombre d'�l�ments
%sinon les solutions des EF sont dans phi et betas
%Les trois premier betas sont les solutions physique. Leur fonction d'onde
%correspondantes sont phi(:,post(1:3)). Chaque colonne de phi correspond a
%une fonction et post te donne la position de la fonction. Donc si tu veux
%la fonction d'onde associ� au deuxi�me beta, ca sera la colonne post(2).
%Sinon, c'est normal si l'amplitude entre les solutions r�elles et d'�F
%sont diff�rentes. Il faut normaliser. Tu divises par le max ou le min la
%solution r�elle et *-1 d�pendament du cas et tu obtiendras la m�me
%expression. c'est normal.

%Sinon, � la toute fin, tu as les valeurs de uEF, duEF et xEF correspondant
%� la solution des �l�ments finis.


% %% Parametres
% m = 2;
% xMin=0;
% xMax=rho*m;
%nEls=m*10; %Doit �tre le m�me multiple que xMax pour qu'un node soit � rho
xn = 0;
L=0;
p = [pDegFM 1];

%% Cr�er le mesh
[nN,nodes,connect,nB,bEls,bPts]=meshEF(xMin,xMax,nEls);

%% degrees of freedom

pDeg=zeros(nEls,1);
pDeg(:,1)=p(1);     %set polynomial degree for each element
pType=zeros(nEls,1);
pType(:,1)=p(2);   %set element type: 1=Lagrangian, 2=hierarchical

[elDof,dFreedom]=dof(nEls,pDeg,connect);


%% Calculer les matrices du probl�me
%Calculer la matrice de gauche d'�l�ments finis
[A] = KUF_A(nEls, xMin, xMax, xn, L, p);
A(end,:)=0.;A(:,end)=0.;A(end,end)=1.;
[As] = A(1:end-1,1:end-1);
%Calculer la matrice de droite d'�l�ments finis
[B] = KUF_B(nEls, xMin, xMax, xn, p);
B(end,:)=0.;B(:,end)=0.;
[Bs] = B(1:end-1,1:end-1);

%% R�soudre le probl�me de valeurs propres
[phi beta2] = eig(As, Bs);
beta = sqrt(beta2);

[sy sx] = size(beta2);
for n = 1:sx
beta2d(n) = beta2(n,n);
end
[beta2s, post] = sort(beta2d,'descend');
betas = sqrt(beta2s);

%% Obtenir un plot et les valeurs de u, du et x
[uEF, duEF, xEF] = postProcFibre(nEls-1,nodes,connect,elDof,dFreedom,pDeg,pType,phi(:,post(VectPropre)), phir, dphir, r, pDegFM, fig);

%% Obtenir la solution exacte
% [r phir betar]=SolutionExacte(400,xMax); %(Numbre de points pour le plot,jusqu'� o� va x)
end
