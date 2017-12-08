function [xiQ,wQ]=gQuad
%sets gaussian quadrature points and weights for orders 1,2,3,4

%order 1
xiQ(1,1)=0;
wQ(1,1)=2;
%order 2
xiQ(1,2)=-1/sqrt(3.);
xiQ(2,2)=-xiQ(1,2);
wQ(1,2)=1;
wQ(2,2)=wQ(1,2);
%order 3
xiQ(1,3)=-sqrt(3./5);
xiQ(2,3)=0;
xiQ(3,3)=-xiQ(1,3);
wQ(1,3)=5./9;
wQ(2,3)=8./9;
wQ(3,3)=wQ(1,3);
%order 4
xiQ(1,4)=-0.861136311594053;
xiQ(2,4)=-0.339981043584856;
xiQ(3,4)=-xiQ(2,4);
xiQ(4,4)=-xiQ(1,4);
wQ(1,4)=0.347854845137454;
wQ(2,4)=0.652145154862546;
wQ(3,4)=wQ(2,4);
wQ(3,4)=wQ(1,4);