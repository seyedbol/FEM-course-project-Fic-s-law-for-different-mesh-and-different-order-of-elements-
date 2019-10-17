function [ say ] = shapefunction( i ) 
global x y ;
if(i==1)
say = (1/4)*(1-x)*(1-y) ;
end
if(i==2)
say = (1/4)*(1+x)*(1-y) ;
end
if(i==3)
say = (1/4)*(1+x)*(1+y) ;
end
if(i==4)
say = (1/4)*(1-x)*(1+y) ;
end

end