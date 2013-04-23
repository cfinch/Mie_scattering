function [s1,s2,qext,qsca,qback,gsca]=mie(x,refrel,angles)
% Calculated based on Mie scattering theory  
% input:
%      x - size parameter =2pi*lambda/radius
%      refrel - refreation index in complext form for example:  1.5+0.02*i;
%      angles - angles at which to compute S1 and S2 function
% output:
%        S1, S2 - funtion which coresponted to phase function
%        Qext - extinction efficiency
%        Qsca - scattering efficiency 
%        Qback -backscatter efficiency
%        gsca- asymmetry parameter

% zatem w sumie jest ich 2*nang-1 bo od 0 do pi

 nmxx=150000;
 nang=length(angles);
 
 s1=zeros(1,nang);     % ilosc katow dla funkcji S1 S2  
 s2=zeros(1,nang);
 d=zeros(1,nmxx);
 pi=zeros(1,nang);
 pi0=zeros(1,nang);
 pi1=ones(1,nang);
 tau=zeros(1,nang);

  pii = 4.*atan(1.);
  dx = x;
  
  drefrl = refrel;
  y = x*drefrl;
  ymod = abs(y);


%    Series expansion terminated after NSTOP terms
%    Logarithmic derivatives calculated from NMX on down

 xstop = x + 4.*x^0.3333 + 2.;
 nmx = max(xstop,ymod) + 15;
 nmx=fix(nmx);
 
% BTD experiment 91/1/15: add one more term to series and compare resu<s
%      NMX=AMAX1(XSTOP,YMOD)+16
% test: compute 7001 wavelen>hs between .0001 and 1000 micron
% for a=1.0micron SiC grain.  When NMX increased by 1, only a single
% computed number changed (out of 4*7001) and it only changed by 1/8387
% conclusion: we are indeed retaining enough terms in series!
      nstop = xstop;
%
      if (nmx > nmxx) %then begin
          'error: nmx > nmxx=', nmxx, ' for |m|x=', ymod
          return
      end
      amu = cos(angles);
      nn = 2*nang - 1;
% Logarithmic derivative D(J) calculated by downward recurrence
% beginning with initial value (0.,0.) at J=NMX
%
      %?d(nmx) = d(0.,0.)
      nn = nmx - 1;
      for n=1: nn   %DO 40 n = 1, nn
          en = nmx - n + 1;
          d(nmx-n) = (en/y) - (1./ (d(nmx-n+1)+en/y));
      end %endfor %40 CONTINUE
%
%*** Riccati-Bessel functions with real argument X
%    calculated by upward recurrence
%
      psi0 = cos(dx);
      psi1 = sin(dx);
      chi0 = -sin(dx);
      chi1 = cos(dx);
      xi1 = psi1-chi1*1i;
      qsca = 0.;
      gsca = 0.;
      p = -1;
      for n=1: nstop  % DO 80 n = 1, nstop
          en = n;
          fn = (2.*en+1.)/ (en* (en+1.));
% for given N, PSI  = psi_n        CHI  = chi_n
%              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
%              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
% Calculate psi_n and chi_n
          psi = (2.*en-1.)*psi1/dx - psi0;
          chi = (2.*en-1.)*chi1/dx - chi0;
          xi = psi-chi*1i;
%
%*** Store previous values of AN and BN for use
%    in computation of g=<cos(theta)>
          if (n > 1) %then begin
              an1 = an;
              bn1 = bn;
          end %endif
%
%*** Compute AN and BN:
          an = (d(n)/drefrl+en/dx)*psi - psi1;
          an = an/ ((d(n)/drefrl+en/dx)*xi-xi1);
          bn = (drefrl*d(n)+en/dx)*psi - psi1;
          bn = bn/ ((drefrl*d(n)+en/dx)*xi-xi1);
%
%*** Augment sums for Qsca and g=<cos(theta)>
          qsca = qsca + (2.*en+1.)* (abs(an)^2+abs(bn)^2);
          gsca = gsca + ((2.*en+1.)/ (en* (en+1.)))* ...
             ( real(an)* real(bn)+imag(an)*imag(bn));

          if (n > 1) %then begin
                     gsca = gsca + ((en-1.)* (en+1.)/en)*...
                    ( real(an1)* real(an)+imag(an1)*imag(an)+...
                     real(bn1)* real(bn)+imag(bn1)*imag(bn));

          end %endif
%
%*** Now calculate scattering intensity pattern
%    First do angles from 0 to 90
          for j=1: nang   %DO 50 j = 1, nang
              jj = 2*nang - j;
              pi(j) = pi1(j);
              tau(j) = en*amu(j)*pi(j) - (en+1.)*pi0(j);
              s1(j) = s1(j) + fn* (an*pi(j)+bn*tau(j));
              s2(j) = s2(j) + fn* (an*tau(j)+bn*pi(j));
          end %endfor % 50     CONTINUE
%
%*** Now do angles greater than 90 using PI and TAU from
%    angles less than 90.
%    P=1 for N=1,3,...% P=-1 for N=2,4,...
%           p = -p;
%           for j=1: nang-1   % DO 60 j = 1, nang - 1
%               jj = 2*nang - j;
%               s1(jj) = s1(jj) + fn*p* (an*pi(j)-bn*tau(j));
%               s2(jj) = s2(jj) + fn*p* (bn*pi(j)-an*tau(j));
%           end %endfor % 60     CONTINUE
          psi0 = psi1;
          psi1 = psi;
          chi0 = chi1;
          chi1 = chi;
          xi1 = psi1-chi1*1i;
%
%*** Compute pi_n for next value of n
%    For each angle J, compute pi_n+1
%    from PI = pi_n , PI0 = pi_n-1
          for j=1: nang   % DO 70 j = 1, nang
              pi1(j) = ((2.*en+1.)*amu(j)*pi(j)- (en+1.)*pi0(j))/...
                      en;
              pi0(j) = pi(j);
           end %endfor %70     CONTINUE
      end %endfor %   80 CONTINUE
%
%*** Have summed sufficient terms.
%    Now compute QSCA,QEXT,QBACK,and GSCA
      gsca = 2.*gsca/qsca;
      qsca = (2./ (dx*dx))*qsca;
      qext = (4./ (dx*dx))* real(s1(1));
      qback = (abs(s1(length(s1))) / dx)^2 / pii;

      ss1=s1;
      ss2=s2;
      clear s1 s2
      a=find(ss1~=0);
      n=max(a);

      s1=ss1(1:n);
      s2=ss2(1:n);
