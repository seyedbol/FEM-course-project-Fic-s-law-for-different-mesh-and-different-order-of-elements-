function [ A ] = intquad(F)
 global x y ;
w=zeros(4,1) ;
p=zeros(4,1) ;
A=0;
w(1,1)=.347854 ; w(2,1)=.652145 ; w(3,1)=.652145 ; w(4,1)=.347854 ;
p(1,1)=-.861136 ; p(2,1)=-.339981 ; p(3,1)=.339981 ; p(4,1)=.861136 ;
for i=1:1:4
    for j=1:1:4
A = A + ( w(i,1)*w(j,1)* subs( F,[x,y],[p(i,1),p(j,1)] ) ) ;

    end
end

end
