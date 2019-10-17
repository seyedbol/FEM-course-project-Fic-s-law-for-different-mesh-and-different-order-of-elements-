function [ C ] = jacobian( X,Y,n ) 

syms x y ;

A1 = (1/4)*(1-x)*(1-y) ;

A2 = (1/4)*(1+x)*(1-y) ;

A3 = (1/4)*(1-x)*(1+y) ;

A4 = (1/4)*(1+x)*(1+y) ;


A = [ diff(A1,x,1) diff(A2,x,1) diff(A3,x,1) diff(A4,x,1) ; diff(A1,y,1) diff(A2,y,1) diff(A3,y,1) diff(A4,y,1) ] ;                              
B = [X(n,1) Y(n,1) ; X(n,3) Y(n,3) ; X(n,7) Y(n,7) ; X(n,9) Y(n,9) ] ;


C= A*B ;





end