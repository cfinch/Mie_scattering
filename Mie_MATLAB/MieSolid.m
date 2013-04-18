(* :Title: MIE Electromagnetic Scattering *)

(* :Context: Mie *)

(* :Summary: Electromagnetic scattering computations for dielectric spheres *)

(* :Copyright: (c) 2003 by Sergio R. Aragon *)

(* :Package Version: 1.0 *)

(* :Sources: Boren & Huffman *)

BeginPackage["MieDefinitions`"]

pi::usage = "pi[n,x,(prec)] is a Legendre type function.
	n = order. Integer >=1.
	x = value of Cos[theta]. Real.
	prec = precision, defaults to 20. "

tau::usage = "tau[n,x,(prec)] is a Legendre type function derivative.
	n = order. Integer >=1.
	x = value of Cos[theta]. Real.
	prec = precision, defaults to 20. "

psi::usage = "psi[n,z,(prec)] is a Ricatti-Bessel function.
	n = order. Integer >=1.
	z = size parameter. Real.
	prec = precision, defaults to 20. "

psip::usage = "psip[n,z] is a Ricatti-Bessel function derivative.
	n = order. Integer >=1.
	z = size parameter. Real.
	 precision is inherited from psi. "

chi::usage = "chi[n,z,(prec)] is a Ricatti-Bessel function.
	n = order. Integer >=1.
	z = size parameter. Real.
	prec = precision, defaults to 20. "

chip::usage = "chip[n,z] is a Ricatti-Bessel function derivative.
	n = order. Integer >=1.
	z = size parameter. Real.
	 precision is inherited from chi. "

zeta::usage = "zeta[n,z] is a Ricatti-Bessel function.
	n = order. Integer >=1.
	z = size parameter. Real.
	 precision is inherited from chi. "

zetap::usage = "zetap[n,z] is a Ricatti-Bessel function derivative.
	n = order. Integer >=1.
	z = size parameter. Real.
	 precision is inherited from chi. "


Begin["`Definitions`"]

(* Legendre Functions *)

pi[n_,x_,Prec_:Automatic]:=
Module[{prec=Prec},	
If[Prec===Automatic, prec = 20];
pi[n,x,prec] = -N[LegendreP[n,1,x]/Sqrt[1-x^2],prec]
]

tau[n_,x_,Prec_:Automatic]:=
Module[{prec=Prec},		
If[Prec===Automatic, prec = 20];
tau[n,x,prec]=-N[(LegendreP[n,2,x]-n(n+1)LegendreP[n,0,x])/2,prec]
]

(* Ricatti-Bessel Functions *)

psi[n_,z_,Prec_:Automatic]:=
Module[{prec=Prec},		
If[Prec===Automatic, prec = 20];
psi[n,z,prec]= N[Sqrt[Pi/2]BesselJ[n+1/2,z]Sqrt[z],prec]
]

psip[n_,z_,Prec_:Automatic]:=
Module[{prec=Prec},		
If[Prec===Automatic, prec = 20];
psip[n,z,prec]= psi[n,z,prec]/z +( n psi[n-1,z,prec]- (n+1) psi[n+1,z,prec])/(2n+1)
]

chi[n_,z_,Prec_:Automatic]:=
Module[{prec=Prec},		
If[Prec===Automatic, prec = 20];
chi[n,z,prec]=-N[Sqrt[Pi/2]BesselY[n+1/2,z]Sqrt[z],prec]
]

chip[n_,z_,Prec_:Automatic]:=
Module[{prec=Prec},		
If[Prec===Automatic, prec = 20];
chip[n,z,prec]=chi[n,z,prec]/z +( n chi[n-1,z,prec]- (n+1) chi[n+1,z,prec])/(2n+1)
]

zeta[n_,z_,Prec_:Automatic]:=
Module[{prec=Prec},		
If[Prec===Automatic, prec = 20];
zeta[n,z,prec]=psi[n,z,prec]-I chi[n,z,prec]
]

zetap[n_,z_,Prec_:Automatic]:=
Module[{prec=Prec},		
If[Prec===Automatic, prec = 20];
zetap[n,z,prec]=psip[n,z,prec]-I chip[n,z,prec]
]

End[]

EndPackage[]

BeginPackage["Mie`","MieDefinitions`"]


a::usage = "a[n,m,z] is the electric Mie scattering coefficient;
	n = the order of the coefficient. n>= 1. Integer;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real."

b::usage = "b[n,m,z] is the magnetic Mie scattering coefficient;
	n = the order the the coefficient. n>= 1. Integer;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real."

S1::usage = "S1[m,z,x] is the vertical scattering amplitude;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real;
	x = value of Cos[theta], theta scattering angle. Real. "

S2::usage = "S2[m,z,x] is the horizontal scattering amplitude;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real;
	x = value of Cos[theta], theta scattering angle. Real."

qscat::usage = "qscat[m,z] is the normalized scattering crossection;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real"

qext::usage = "qext[m,z_] is the total normalized extinction crossection;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real."

qabs::usage = "qabs[m,z] is the normalized absorption crossection;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real."

albedo::usage = "albedo[m,z] is the fraction of extinction due to scattering;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real."

asymmetry::usage = "asymmetrys[m,z] is <Cos[theta]>, theta = scattering angle;
	m = the relative index of refraction of the sphere (can be complex);
	z = kr, where r is the sphere radius, k = 2Pi/lambda.  Real."

$MaxExtraPrecision = 1000;

Begin["Mie`SolidDielectric`"]

prec=20;
a[n_,m_,z_]:=(m psi[n,m z,prec] psip[n,z,prec]-psi[n,z,prec] psip[n,m z,prec])/(m psi[n,m z,prec] zetap[n,z,prec]-zeta[n,z,prec] psip[n,m z,prec]);
b[n_,m_,z_]:=(psi[n,m z,prec] psip[n,z]-m psi[n,z,prec] psip[n,m z])/(psi[n,m z,prec] zetap[n,z]-m zeta[n,z] psip[n,m z]);

S1[m_,z_,x_] := Sum[((2n + 1)/(n(n + 1))(a[n, m, z]pi[n, x,prec] + b[n, m, z]tau[n, x,prec])),{n,1,IntegerPart[z + 4 z^(1/3) + 2]}];

S2[m_,z_,x_] := Sum[((2n + 1)/(n(n + 1))(b[n, m, z]pi[n, x,prec] + a[n, m, z]tau[n, x,prec])),{n,1,IntegerPart[z + 4 z^(1/3) + 2]}];

qscat[m_,z_] := 
2/(z^2)Sum[ ((2n+1)(Abs[a[n,m,z]]^2+Abs[b[n,m,z]]^2)),{n,1,IntegerPart[z + 4 z^(1/3) + 2]}];

qext[m_,z_] := 
2/(z^2)Sum[ ((2n + 1)Re[a[n, m, z] + b[n, m, z]]),{n,1,IntegerPart[z + 4 z^(1/3) + 2]}];

qabs[m_,z_] := qext[m, z] - qscat[m, z];

albedo[m_,z_] := qscat[m,z]/qext[m,z];

asymmetry[m_,z_] := (4/(z^2 qscat[z,m])) Sum[(n(n+2)( Re[a[n,m,z] Conjugate[a[n+1,m,z]]
	+b[n,m,z] Conjugate[b[n+1,m,z]]])+(2n+1)Re[a[n,m,z]Conjugate[b[n,m,z]]]/n)/(n+1),{n,1,IntegerPart[z + 4 z^(1/3) + 2]}];

End[]

EndPackage[]