function result = Mie(m, x)
% Computation of Mie Efficiencies for given 
% complex refractive-index ratio m=m'+im" 
% and size parameter x=k0*a, where k0= wave number in ambient 
% medium, a=sphere radius, using complex Mie Coefficients
% an and bn for n=1 to nmax,
% s. Bohren and Huffman (1983) BEWI:TDD122, p. 103,119-122,477.
% Result: m', m", x, efficiencies for extinction (qext), 
% scattering (qsca), absorption (qabs), backscattering (qb), 
% asymmetry parameter (asy=<costeta>) and (qratio=qb/qsca).
% Uses the function "Mie_ab" for an and bn, for n=1 to nmax.
% C. Mätzler, May 2002, revised July 2002.

if x==0                 % To avoid a singularity at x=0
    result=[0 0 0 0 0 1.5];
elseif x>0              % This is the normal situation
    nmax=round(2+x+4*x.^(1/3));
    n1=nmax-1;
    n=(1:nmax);cn=2*n+1; c1n=n.*(n+2)./(n+1); c2n=cn./n./(n+1);
    x2=x.*x;
    f=mie_ab(m,x);
    anp=(real(f(1,:))); anpp=(imag(f(1,:)));
    bnp=(real(f(2,:))); bnpp=(imag(f(2,:)));
    g1(1:4,nmax)=[0; 0; 0; 0]; % displaced numbers used for
    g1(1,1:n1)=anp(2:nmax);    % asymmetry parameter, p. 120
    g1(2,1:n1)=anpp(2:nmax);
    g1(3,1:n1)=bnp(2:nmax);
    g1(4,1:n1)=bnpp(2:nmax);   
    dn=cn.*(anp+bnp);
    q=sum(dn);
    qext=2*q/x2;
    en=cn.*(anp.*anp+anpp.*anpp+bnp.*bnp+bnpp.*bnpp);
    q=sum(en);
    qsca=2*q/x2;
    qabs=qext-qsca;
    fn=(f(1,:)-f(2,:)).*cn;
    gn=(-1).^n;
    f(3,:)=fn.*gn;
    q=sum(f(3,:));
    qb=q*q'/x2;
    asy1=c1n.*(anp.*g1(1,:)+anpp.*g1(2,:)+bnp.*g1(3,:)+bnpp.*g1(4,:));
    asy2=c2n.*(anp.*bnp+anpp.*bnpp);
    asy=4/x2*sum(asy1+asy2)/qsca;
    qratio=qb/qsca;
    result=[qext qsca qabs qb asy qratio];
end;