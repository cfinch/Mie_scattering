function result=Mie_pt(u,nmax)
% pi_n and tau_n, -1 <= u <= 1, n1 integer from 1 to nmax 
% angular functions used in Mie Theory
% Bohren and Huffman (1983), p. 94 - 95

p(1)=1; 
t(1)=u;
p(2)=3*u; 
t(2)=3*cos(2*acos(u));
for n1=3:nmax,
    p1=(2*n1-1)./(n1-1).*p(n1-1).*u;
    p2=n1./(n1-1).*p(n1-2);
    p(n1)=p1-p2;
    t1=n1*u.*p(n1);
    t2=(n1+1).*p(n1-1);
    t(n1)=t1-t2;
end;

result=[p;t];
