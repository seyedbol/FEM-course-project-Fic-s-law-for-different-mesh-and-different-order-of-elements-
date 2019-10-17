function [ a ] = K( x,y )
a = (10^-9)+( (4*10^-9) * sin((pi*x)/2) * y) ;
end