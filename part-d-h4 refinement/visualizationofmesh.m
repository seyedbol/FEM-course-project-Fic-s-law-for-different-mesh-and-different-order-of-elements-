clc
clear all
close all

load Elements.txt
load Nodes.txt

nrelm=size(Elements,1);
nnod=size(Nodes,1);

Edof=zeros(nrelm,9);
Ex=zeros(nrelm,4); 
Ey=Ex;
Elm(:,1)=sort(1:nrelm);

for i=1:nrelm
 
    Elm(i,2:5)=Elements(i,6:9);  
    Ex(i,:)=[Nodes(Elm(i,2),2) Nodes(Elm(i,3),2) Nodes(Elm(i,4),2) Nodes(Elm(i,5),2)];
    Ey(i,:)=[Nodes(Elm(i,2),3) Nodes(Elm(i,3),3) Nodes(Elm(i,4),3) Nodes(Elm(i,5),3)];
    
end




figure(1);
hold on;
 for i=1:nrelm     %%%length(Edof(:,1))
    plot(Ex(i,[1:end 1]),Ey(i,[1:end 1]));
end
axis equal