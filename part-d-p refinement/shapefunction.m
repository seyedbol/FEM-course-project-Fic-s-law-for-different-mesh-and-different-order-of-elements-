function [ say ] = shapefunction( i ) 
global x y ;

if(i==1)
say = (1/4)*( (x^2)-x )*( (y^2)-y ) ;
end
if(i==2)
say = (1/2)*( 1-(x^2))*( (y^2)-y ) ;
end
if(i==3)
say = (1/4)*( (x^2)+x )*( (y^2)-y )  ;
end
if(i==4)
say = (1/2)*( (x^2)-x )*( 1-(y^2) ) ;
end
if(i==5)
say = ( 1-(x^2))*( 1-(y^2) ) ;
end
if(i==6)
say = (1/2)*( (x^2)+x )*( 1-(y^2) ) ;
end
if(i==7)
say = (1/4)*( (x^2)-x )*( (y^2)+y ) ;
end
if(i==8)
say = (1/2)*( 1-(x^2))*( (y^2)+y  ) ;
end
if(i==9)
say = (1/4)*( (x^2)+x )*( (y^2)+y ) ;
end

end