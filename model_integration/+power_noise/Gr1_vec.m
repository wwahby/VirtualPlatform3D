% Green function
function Gr1=Gr1_vec(x,y,w,v,lambda,a)
%global testy

M=1000;

jj = 0:M-1;
p = pi*jj/a;
q = pi*jj/a;
beta = sqrt(p.^2-lambda);
ep = 2*ones(1,M);
ep(1) = 1;

% Decide whether we use y or v in the temp accumulator
if (y>=v)
    yy = y;
    vv = v;
else
    yy = v;
    vv = y;
end

temp = zeros(1,M);
ind1 = 101;
ep1 = ep(1:ind1);
b1 = beta(1:ind1);
p1 = p(1:ind1);

ep2 = ep(ind1+1:end);
b2 = beta(ind1+1:end);
p2 = p(ind1+1:end);


temp(1:ind1) = 1./a.*(ep1.*cos(p1*x).*cos(p1*w)).*cosh(b1*vv).*cosh(b1*(a-yy)) ./ ( b1.*sinh(b1*a))  ;
temp(ind1+1:end) = 1./a.*(ep2.*cos(p2*x).*cos(p2*w)).*exp( b2*vv + b2*(a-yy) - b2*a) ./ ( 2*b2 )  ;


% if (y>=v)
%     temp(1:101) = 1/a.*(ep.*cos(p*x).*cos(p*w)) / ( beta.*sinh(beta*a).*cosh(beta*v).*cosh(beta*(a-y)) );
%     temp(102:end) = 1/a.*(ep.*cos(p*x).*cos(p*w)) / ( 2*beta.*exp( beta*v + beta*(a-y) - beta*a) );
% else
%     temp(1:101) = 1/a.*(ep.*cos(p*x).*cos(p*w)) / ( beta.*sinh(beta*a).*cosh(beta*y).*cosh(beta*(a-v)) );
%     temp(102:end) = 1/a.*(ep.*cos(p*x).*cos(p*w)) / ( 2*beta.*exp(beta*y + beta*(a-v) - beta*a) );
% end

Gr1 = sum(temp);
%testy = abs(cumsum(temp));

% 
% for j=0:M-1
%     p(j+1)=pi*j/a;
%     q(j+1)=pi*j/a;
%     beta(j+1)=sqrt(p(j+1)*p(j+1)-lambda);
% end
% 
% ep=zeros(1,M);
% ep(1)=1;
% for j=1:M-1
%     ep(j+1)=2;
% end
% 
% testy=zeros(1,M);
% testx=zeros(1,M);
% 
% temp=0;
% testv=0;
% testn=1;

% 
% for n=0:M-1
%     if y>=v
%         if n<=100
%         temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)/sinh(beta(n+1)*a)*cosh(beta(n+1)*v)*cosh(beta(n+1)*(a-y));
%         end
%         
%         if n>100
%             temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)*exp(beta(n+1)*v+beta(n+1)*(a-y)-beta(n+1)*a)/2;
%         end
% 
%     end
%     if y<v
%         if n<=100
%           temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)/sinh(beta(n+1)*a)*cosh(beta(n+1)*y)*cosh(beta(n+1)*(a-v));
%         end
%         if n>100
%           temp=temp+1/a*(ep(n+1)*cos(p(n+1)*x)*cos(p(n+1)*w))/beta(n+1)*exp(beta(n+1)*y+beta(n+1)*(a-v)-beta(n+1)*a)/2; 
%        end
%     end
%      
%     testy(testn)=abs(temp);
%     testx(testn)=testn;
%     testn=testn+1;
%     
% end
% 
% Gr1=temp;
% end
