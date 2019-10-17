function [ F2 ] = mapx( X,n ) 
global x y ;


A1 = (1/4)*(1-x)*(1-y) ;

A2 = (1/4)*(1+x)*(1-y) ;

A3 = (1/4)*(1-x)*(1+y) ;

A4 = (1/4)*(1+x)*(1+y) ;


F2 = (A1*X(n,1))+(A2*X(n,3))+(A3*X(n,7))+(A4*X(n,9)) ;

end


