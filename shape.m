function [N,dN]=shape(xi,ne,pDeg,pType)
%calculates the values of shape functions N and their derivatives
%dN/d(xi) at specified values of xi

%can be linear or quadratic, Lagrangian or hierarchical for each element

if pDeg(ne)==1  %linear (only one form)
    N(1)=.5*(1-xi);
    N(2)=.5*(1+xi);
    dN(1)=-.5;
    dN(2)=.5;
end
if pDeg(ne)==2;  %quadratic
    if pType(ne)==1;   %Lagrangian
        N(1)=.5*xi*(xi-1.);
        N(3)=1-xi^2;
        N(2)=.5*xi*(xi+1.);
        dN(1)=xi-.5;
        dN(3)=-2*xi;
        dN(2)=xi+.5;
    end
    if pType(ne)==2;   %hierarchical
        N(1)=.5*(1-xi);
        N(2)=.5*(1+xi);
        N(3)=1-xi^2;
        dN(1)=-.5;
        dN(2)=.5;
        dN(3)=-2*xi;
    end
end

end

    