% Green function
function Gr1=Gr1(x,y,w,v,lamda,a)
global testy

M=1000;
p=zeros(1,M);
q=zeros(1,M);
beta=zeros(1,M);

for j=0:M-1
    p(j+1)=pi*j/a;
    q(j+1)=pi*j/a;
    beta(j+1)=sqrt(p(j+1)*p(j+1)-lamda);
end

ep=zeros(1,M);
ep(1)=1;
for j=1:M-1
    ep(j+1)=2;
end

testy=zeros(1,M);
testx=zeros(1,M);

temp=0;
testv=0;
testn=1;
for n=0:M-1
    if y>=v
        if n<=100
        temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)/sinh(beta(n+1)*a)*cosh(beta(n+1)*v)*cosh(beta(n+1)*(a-y));
        end
        
        if n>100
            temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)*exp(beta(n+1)*v+beta(n+1)*(a-y)-beta(n+1)*a)/2;
        end
        %sinh(beta(n+1)*a);
        %n
        %beta(n+1)
        %exp(beta(n+1)*a)
        %exp(-beta(n+1)*a)
        %cosh(beta(n+1)*v)
    end
    if y<v
        if n<=100
          temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)/sinh(beta(n+1)*a)*cosh(beta(n+1)*y)*cosh(beta(n+1)*(a-v));
        end
        if n>100
          temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)*exp(beta(n+1)*y+beta(n+1)*(a-v)-beta(n+1)*a)/2; 
       end
     end
        testy(testn)=abs(temp);
        testx(testn)=testn;
        testn=testn+1;
    end
%end
%plot(testx,testy)
Gr1=temp;
end
