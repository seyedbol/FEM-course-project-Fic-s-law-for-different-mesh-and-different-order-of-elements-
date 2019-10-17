clc ;
clear all ;


%%%%%%%%%%%%%coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=-9.81 ; miu=0.001 ; ro=1000 ; p0=121300 ; pinfinity=101300 ; 
global x y Ex Ey connectivity connectivitynew ;
syms x y ;
global nnod nrelm ;

%%%%%%%%%%%organizing elements and nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Elements.txt
load Nodes.txt

nrelm=size(Elements,1);
nnod=size(Nodes,1);
Ex=zeros(nrelm,9);
Ey=zeros(nrelm,9);
connectivity(:,1)=sort(1:nrelm);

for i=1:nrelm
 
    connectivity(i,1:9)=Elements(i,6:14); %%%%creating connectivity matrix 
    Ex(i,:)=[Nodes(connectivity(i,1),2) Nodes(connectivity(i,2),2) Nodes(connectivity(i,3),2) Nodes(connectivity(i,4),2) Nodes(connectivity(i,5),2) Nodes(connectivity(i,6),2) Nodes(connectivity(i,7),2) Nodes(connectivity(i,8),2) Nodes(connectivity(i,9),2)];
    Ey(i,:)=[Nodes(connectivity(i,1),3) Nodes(connectivity(i,2),3) Nodes(connectivity(i,3),3) Nodes(connectivity(i,4),3) Nodes(connectivity(i,5),3) Nodes(connectivity(i,6),3) Nodes(connectivity(i,7),3) Nodes(connectivity(i,8),3) Nodes(connectivity(i,9),3)];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%new connectivity

connectivitynew=zeros(nrelm,9) ;
connectivitynew(:,1)=connectivity(:,2) ;
connectivitynew(:,2)=connectivity(:,6) ;
connectivitynew(:,3)=connectivity(:,3) ;
connectivitynew(:,4)=connectivity(:,5) ;
connectivitynew(:,5)=connectivity(:,9) ;
connectivitynew(:,6)=connectivity(:,7) ;
connectivitynew(:,7)=connectivity(:,1) ;
connectivitynew(:,8)=connectivity(:,8) ;
connectivitynew(:,9)=connectivity(:,4) ;
for i=1:nrelm
 
    Ex(i,:)=[Nodes(connectivitynew(i,1),2) Nodes(connectivitynew(i,2),2) Nodes(connectivitynew(i,3),2) Nodes(connectivitynew(i,4),2) Nodes(connectivitynew(i,5),2) Nodes(connectivitynew(i,6),2) Nodes(connectivitynew(i,7),2) Nodes(connectivitynew(i,8),2) Nodes(connectivitynew(i,9),2)];
    Ey(i,:)=[Nodes(connectivitynew(i,1),3) Nodes(connectivitynew(i,2),3) Nodes(connectivitynew(i,3),3) Nodes(connectivitynew(i,4),3) Nodes(connectivitynew(i,5),3) Nodes(connectivitynew(i,6),3) Nodes(connectivitynew(i,7),3) Nodes(connectivitynew(i,8),3) Nodes(connectivitynew(i,9),3)];
    
end
%%%%%%%%%%%%%%%%%creating local matrix (linear 4*4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localmatrix_l = zeros(9,9) ;%%left
localmatrix_r = zeros(9,1) ;%%right
globalmatrix_l = sparse(nnod,nnod) ;%%left
globalmatrix_r = zeros(nnod,1) ;%%right


 for n=1:1:nrelm%%%%swipe the elements
     
     for i=1:1:9%%%swipe the rows of local matrix
  
  sayi=shapefunction(i)  ;%%shape function
  Hi=inv(jacobian(Ex,Ey,n)) ;
  F2=(  K( mapx( Ex,n), mapy( Ey,n) ) / miu ) * g * ro * ( Hi(2,1)*diff(sayi,x,1)+Hi(2,2)*diff(sayi,y,1) )*det(jacobian(Ex,Ey,n));
  %calculating the integral with the help of quadrature points(RIGHT HAND SIDE)
  localmatrix_r(i,1) =intquad(F2) ;
  %%putting local matrix of right hand side into global matrix
  globalmatrix_r(connectivitynew(n,i)) = globalmatrix_r(connectivitynew(n,i)) + localmatrix_r(i,1) ;

      
        for j=1:1:9%%%swipe the columns of local matrix
              
  %%%%%calculating shape functions&&calculating F(zeta) ( LEFT HAND SIDE )
  sayj=shapefunction(j)  ;%%%shape function for pressure
  F1=( K( mapx( Ex,n), mapy( Ey,n) ) / miu )*det(jacobian(Ex,Ey,n))...
  * ( ( ( Hi(1,1)*diff(sayi,x,1) )+( Hi(1,2)*diff(sayi,y,1) ) )*( ( Hi(1,1)*diff(sayj,x,1) )+( Hi(1,2)*diff(sayj,y,1) ) )+...
  ( ( Hi(2,1)*diff(sayi,x,1) )+( Hi(2,2)*diff(sayi,y,1) ) )*( ( Hi(2,1)*diff(sayj,x,1) )+( Hi(2,2)*diff(sayj,y,1) )  )  );            
 %calculating the integral with the help of quadrature points(LEFT HAND SIDE) 
  localmatrix_l(i,j) = intquad(F1) ; 
 %%putting local matrix of right left side into global matrix
  globalmatrix_l(connectivitynew(n,i),connectivitynew(n,j)) = globalmatrix_l(connectivitynew(n,i),connectivitynew(n,j)) + localmatrix_l(i,j) ;
 
  %%%%%%%%%%%%%%
        end 
     end
     
 end  
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Drichlet boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:nnod     
    %%%inlet
    
   if(i==7||i==6)
      globalmatrix_l(i,:)=0 ;
     globalmatrix_l(i,i)=1 ;
      globalmatrix_r(i,1)=p0 ;
   end
  if(i>7&&i<13)
          globalmatrix_l(i,:)=0 ;  
          globalmatrix_l(i,i)=1 ;
          globalmatrix_r(i,1)=p0 ;     
  end
  
  %%%outlet
    if(i==4||i==3)
      globalmatrix_l(i,:)=0 ;
      globalmatrix_l(i,i)=1 ;
      globalmatrix_r(i,1)=pinfinity ;  
    end
   if(i>38&&i<44)
          globalmatrix_l(i,:)=0 ;  
          globalmatrix_l(i,i)=1 ;
          globalmatrix_r(i,1)=pinfinity ;
   end 


 end

%%%%%%%%%%%%%%%%%%%%%%%post processing
global Pressure ;
global profile ;   
Pressure = globalmatrix_l\globalmatrix_r ;
profile = zeros(nrelm,9) ;
   
    for iel=1:nrelm  
        nd=connectivitynew(iel,1:end);         % extract connected node for (iel)-th element
        profile(iel,:) = (Pressure(nd)) ;         % extract component value of the node 
    end
    figure(1)
    for i=1:nrelm    
    fill(Ex(i,[1 7 9 3 1]),Ey(i,[1 7 9 3  1]),profile(i,[1 7 9 3 1]));
    hold on ;   
    end
    colormap 'jet'
    title('pressure')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
VX=zeros(nnod,1) ;
VY=zeros(nnod,1) ;
Vfinal=zeros(nnod,1) ;
%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:nnod
    counter=0 ;
    Vx=0;
    Vy=0;
    for n=1:1:nrelm
        for j=1:1:9
if( connectivitynew(n,j)== k)
counter=counter+1 ;
A = (profile(n,1)*shapefunction(1))+(profile(n,2)*shapefunction(2))+(profile(n,3)*shapefunction(3))+(profile(n,4)*shapefunction(4))...
    +( profile(n,5)*shapefunction(5))+(profile(n,6)*shapefunction(6))+(profile(n,7)*shapefunction(7))+(profile(n,8)*shapefunction(8))...
    +(profile(n,9)*shapefunction(9)) ;
    
Hi=inv(jacobian(Ex,Ey,n)) ;
Xvelocity= -(  K( Ex(n,j),Ey(n,j))/  miu   )* ( Hi(1,1)*diff(A,x,1)+Hi(1,2)*diff(A,y,1) ) ;   
Yvelocity= -(  K( Ex(n,j),Ey(n,j)) / miu )* ( (Hi(2,1)*diff(A,x,1)+Hi(2,2)*diff(A,y,1)) - (ro*g) ) ; 
%%%%%%%%%%%%%%
if(j==1)
Vx=Vx+subs(Xvelocity,[x,y],[-1,-1]);
Vy=Vy+subs(Yvelocity,[x,y],[-1,-1]);
end
if(j==2)
Vx=Vx+subs(Xvelocity,[x,y],[0,-1]);
Vy=Vy+subs(Yvelocity,[x,y],[0,-1]);
end
if(j==3)
Vx=Vx+subs(Xvelocity,[x,y],[1,-1]);
Vy=Vy+subs(Yvelocity,[x,y],[1,-1]);
end
if(j==4)
Vx=Vx+subs(Xvelocity,[x,y],[-1,0]);
Vy=Vy+subs(Yvelocity,[x,y],[-1,0]);
end
if(j==5)
Vx=Vx+subs(Xvelocity,[x,y],[0,0]);
Vy=Vy+subs(Yvelocity,[x,y],[0,0]);
end
if(j==6)
Vx=Vx+subs(Xvelocity,[x,y],[1,0]);
Vy=Vy+subs(Yvelocity,[x,y],[1,0]);
end
if(j==7)
Vx=Vx+subs(Xvelocity,[x,y],[-1,1]);
Vy=Vy+subs(Yvelocity,[x,y],[-1,1]);
end
if(j==8)
Vx=Vx+subs(Xvelocity,[x,y],[0,1]);
Vy=Vy+subs(Yvelocity,[x,y],[0,1]);
end
if(j==9)
Vx=Vx+subs(Xvelocity,[x,y],[1,1]);
Vy=Vy+subs(Yvelocity,[x,y],[1,1]);
end
%%%%%%%

 end           
        end

    end
 VX(k,1)=Vx/counter  ;
 VY(k,1)=Vy/counter  ;
 Vfinal(k,1)=sqrt( ( VX(k,1)^2 )+(VY(k,1)^2) ) ;
end

%%%%%%%%%%%%%%%%%post processing
profile2 = zeros(nrelm,9) ;
profile3 = zeros(nrelm,9) ;   
    for iel=1:nrelm  
        nd=connectivitynew(iel,1:end);         % extract connected node for (iel)-th element
        profile2(iel,:) = ((VY(nd))) ;
        profile3(iel,:) = ((VX(nd))) ;  
        % extract component value of the node 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%   
figure (2)    
for i=1:nrelm    
    fill(Ex(i,[1 7 9 3 1]),Ey(i,[1 7 9 3  1]),profile2(i,[1 7 9 3 1]));
    hold on ;   
end
    colormap 'jet'
      title('VY')
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      figure(3)
      for i=1:nrelm    
    fill(Ex(i,[1 7 9 3 1]),Ey(i,[1 7 9 3  1]),profile3(i,[1 7 9 3 1]));
    hold on ;   
      end
    colormap 'jet'
    title('VX')
      