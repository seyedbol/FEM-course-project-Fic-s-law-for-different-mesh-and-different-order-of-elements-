%%%%%%%%%load the data that I saved from main code
load Pressure.mat 
load profile.mat
load Ex.mat 
load Ey.mat
load connectivity 
load Nodes.txt
load Elements.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%coefficients
g=-9.81 ; miu=0.001 ; ro=1000 ; p0=121300 ; u0=-.001 ; pinfinity=101300 ; hmaster=2 ;
global x y ;
syms x y ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%defining matrixes
nnod=size(Nodes,1);
nrelm=size(Elements,1);
VX=zeros(nnod,1) ;
VY=zeros(nnod,1) ;
Vfinal=zeros(nnod,1) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
profile2 = zeros(nrelm,4) ;
   
    for iel=1:nrelm  
        nd=connectivity(iel,1:end);         % extract connected node for (iel)-th element
        profile2(iel,:) = ((VY(nd))) ;         % extract component value of the node 
    end
    
    
    for i=1:nrelm    
    fill(Ex(i,[1:end 1]),Ey(i,[1:end 1]),profile2(i,[1:end 1]));
    hold on ;   
    end
    colormap 'jet'
    
    
   
    
    
    
    