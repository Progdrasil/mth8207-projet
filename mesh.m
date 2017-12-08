function [nN,nodes,connect,nB,bEls,bPts]=mesh(xMin,xMax,nEls)
%inputs domain boundary points, # of elements; outputs # of nodes,
%their coordinates, connectivity matrix, # of boundary points

nN=nEls+1;  %# of nodes
nodes=(xMin: (xMax-xMin)/nEls: xMax)'; %coordinates of element edges
connect=zeros(nEls,2);  %init connectivity matrix
i=1;
for j=1:nEls
    connect(j,1)=i;
    i=i+1;
    connect(j,2)=i;
end

%boundary conditions
nB=2;               %# of boundary elements
bEls=zeros(2,1);    %boundary elements
bEls(1)=1;
bEls(2)=nEls;
bPts=zeros(2,1);  %boundary point nodes
bPts(1)=1;
bPts(2)=2;
end