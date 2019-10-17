clc ;
clear all ;


%%%%%%%%%%%%%coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=-9.81 ; miu=0.001 ; ro=1000 ; p0=121300 ; pinfinity=101300 ; 
global x y Ex Ey connectivity ;
syms x y ;
global nnod nrelm ;

%%%%%%%%%%%organizing elements and nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Elements.txt
load Nodes.txt

nrelm=size(Elements,1);
nnod=size(Nodes,1);
Ex=zeros(nrelm,4);
Ey=zeros(nrelm,4);
connectivity(:,1)=sort(1:nrelm);

for i=1:nrelm
 
    connectivity(i,1:4)=Elements(i,6:9); %%%%creating connectivity matrix 
    Ex(i,:)=[Nodes(connectivity(i,1),2) Nodes(connectivity(i,2),2) Nodes(connectivity(i,3),2) Nodes(connectivity(i,4),2)];
    Ey(i,:)=[Nodes(connectivity(i,1),3) Nodes(connectivity(i,2),3) Nodes(connectivity(i,3),3) Nodes(connectivity(i,4),3)];
    
end



%%%%%%%%%%%%%%%%%creating local matrix (linear 4*4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localmatrix_l = zeros(4,4) ;%%left
localmatrix_r = zeros(4,1) ;%%right
globalmatrix_l = zeros(nnod,nnod) ;%%left
globalmatrix_r = zeros(nnod,1) ;%%right


 for n=1:1:nrelm%%%%swipe the elements
     
     for i=1:1:4%%%swipe the rows of local matrix
  
  sayi=shapefunction(i)  ;%%shape function
  Hi=inv(jacobian(Ex,Ey,n)) ;
  F2=(  K( mapx( Ex,n), mapy( Ey,n) ) / miu ) * g * ro * ( Hi(2,1)*diff(sayi,x,1)+Hi(2,2)*diff(sayi,y,1) )*det(jacobian(Ex,Ey,n));
  %calculating the integral with the help of quadrature points(RIGHT HAND SIDE)
  localmatrix_r(i,1,n) =intquad(F2) ;
  %%putting local matrix of right hand side into global matrix
  globalmatrix_r(connectivity(n,i)) = globalmatrix_r(connectivity(n,i)) + localmatrix_r(i,1) ;

      
        for j=1:1:4%%%swipe the columns of local matrix
              
  %%%%%calculating shape functions&&calculating F(zeta) ( LEFT HAND SIDE )
  sayj=shapefunction(j)  ;%%%shape function for pressure
  F1=( K( mapx( Ex,n), mapy( Ey,n) ) / miu )*det(jacobian(Ex,Ey,n))...
  * ( ( ( Hi(1,1)*diff(sayi,x,1) )+( Hi(1,2)*diff(sayi,y,1) ) )*( ( Hi(1,1)*diff(sayj,x,1) )+( Hi(1,2)*diff(sayj,y,1) ) )+...
  ( ( Hi(2,1)*diff(sayi,x,1) )+( Hi(2,2)*diff(sayi,y,1) ) )*( ( Hi(2,1)*diff(sayj,x,1) )+( Hi(2,2)*diff(sayj,y,1) )  )  );            
 %calculating the integral with the help of quadrature points(LEFT HAND SIDE) 
  localmatrix_l(i,j) = intquad(F1) ; 
 %%putting local matrix of right left side into global matrix
  globalmatrix_l(connectivity(n,i),connectivity(n,j)) = globalmatrix_l(connectivity(n,i),connectivity(n,j)) + localmatrix_l(i,j) ;
 
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

%%%%%%%%%%%%%%%%%%%%%%%%post processing the pressure 
global Pressure ;
global profile ;   
Pressure = globalmatrix_l\globalmatrix_r ;
profile = zeros(nrelm,4) ;
   
    for iel=1:nrelm  
        nd=connectivity(iel,1:end);         % extract connected node for (iel)-th element
        profile(iel,:) = (Pressure(nd)) ;         % extract component value of the node 
    end
    figure(1)
    for i=1:nrelm    
    fill(Ex(i,[1:end 1]),Ey(i,[1:end 1]),profile(i,[1:end 1]));
    hold on ;   
    end
    colormap 'jet'
    title('pressure')
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%velocity calculation
VX=zeros(nnod,1) ;
VY=zeros(nnod,1) ;
Vfinal=zeros(nnod,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:nnod
    counter=0 ;
    Vx=0;
    Vy=0;
    for n=1:1:nrelm
        for j=1:1:4
if( connectivity(n,j)== k )
counter=counter+1 ;
A = (profile(n,1)*shapefunction(1))+(profile(n,2)*shapefunction(2))+(profile(n,3)*shapefunction(3))+(profile(n,4)*shapefunction(4)) ;
Hi=inv(jacobian(Ex,Ey,n)) ;
Xvelocity= -(  K( Ex(n,j),Ey(n,j))/miu)* ( Hi(1,1)*diff(A,x,1)+Hi(1,2)*diff(A,y,1) ) ;   
Yvelocity= -(  K( Ex(n,j),Ey(n,j)) / miu )* ( (Hi(2,1)*diff(A,x,1)+Hi(2,2)*diff(A,y,1)) - (ro*g) ) ; 
if(j==1)
Vx=Vx+subs(Xvelocity,[x,y],[-1,-1]);
Vy=Vy+subs(Yvelocity,[x,y],[-1,-1]);
end
if(j==2)
Vx=Vx+subs(Xvelocity,[x,y],[1,-1]);
Vy=Vy+subs(Yvelocity,[x,y],[1,-1]);
end
if(j==3)
Vx=Vx+subs(Xvelocity,[x,y],[1,1]);
Vy=Vy+subs(Yvelocity,[x,y],[1,1]);
end
if(j==4)
Vx=Vx+subs(Xvelocity,[x,y],[-1,1]);
Vy=Vy+subs(Yvelocity,[x,y],[-1,1]);
end
 end           
        end

    end
 VX(k,1)=Vx/counter  ;
 VY(k,1)=Vy/counter  ;
 Vfinal(k,1)=sqrt( ( VX(k,1)^2 )+(VY(k,1)^2) ) ;
end



%%%%%%%%%%%%%%%%%post processing

figure(2)
profile2 = zeros(nrelm,4) ;
profile3 = zeros(nrelm,4) ;   
    for iel=1:nrelm  
        nd=connectivity(iel,1:end);         % extract connected node for (iel)-th element
        profile2(iel,:) = ((VY(nd))) ;
        % extract component value of the node 
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%   
    for i=1:nrelm    
    fill(Ex(i,[1:end 1]),Ey(i,[1:end 1]),profile2(i,[1:end 1]));
    hold on ;   
    end
    colormap 'jet'
    title('Vy')
    %%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%
       for iel=1:nrelm  
        nd=connectivity(iel,1:end);         % extract connected node for (iel)-th element
        profile3(iel,:) = ((VX(nd))) ;         % extract component value of the node 
       end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
        for i=1:nrelm    
    fill(Ex(i,[1:end 1]),Ey(i,[1:end 1]),profile3(i,[1:end 1]));
    hold on ;   
    end
    colormap 'jet'
    title('Vx')
