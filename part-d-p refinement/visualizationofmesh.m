clc
clear all
close all

load Elements.txt
load Nodes.txt

nrelm=size(Elements,1);
nnod=size(Nodes,1);

Ex=zeros(nrelm,9); 
Ey=Ex;
Elm(:,1)=sort(1:nrelm);

for i=1:nrelm
 
    Elm(i,2:10)=Elements(i,6:14);  
    Ex(i,:)=[Nodes(Elm(i,2),2) Nodes(Elm(i,3),2) Nodes(Elm(i,4),2) Nodes(Elm(i,5),2) Nodes(Elm(i,6),2) Nodes(Elm(i,7),2) Nodes(Elm(i,8),2) Nodes(Elm(i,9),2) Nodes(Elm(i,10),2)];
    Ey(i,:)=[Nodes(Elm(i,2),3) Nodes(Elm(i,3),3) Nodes(Elm(i,4),3) Nodes(Elm(i,5),3) Nodes(Elm(i,6),3) Nodes(Elm(i,7),3) Nodes(Elm(i,8),3) Nodes(Elm(i,9),3) Nodes(Elm(i,10),3)];
    
end




figure(1);
hold on;
 for i=1:nrelm     %%%length(Edof(:,1))
    plot(Ex(i,[1:4 1]),Ey(i,[1:4 1]));
end
axis equal