var
np ng u m Am v vp vg mu fp
fg chi nu n tauc lamb psi I r Y
K N z Op Og wp tauw wg J pm
y GK wpstar pistar pi F1 F2 b gc k
gi i T C taur gs Ag MUg Af Uc 
thet hp hg ;

parameters
kapp chig gamm zet d bet delt alph lp lg et
omegaw omegap sigmaa sigmab xi ha hb
rhogc rhogi rhoz rhotauw rhotauc rhotaur 
phi rhoMUg rhovg rhoi phipi phiy Ga rhoAf rhoAm rhoAg
rhonu rhochi gi_Y gc_Y b_Y T_Y z_wp
;
% basic param
delt=0.025;
bet=0.97;
d=0.6;
Ga=0.05;
% price param
omegap=0.5;
omegaw=0.5;
et=11;
%param wage and match
kapp=0.5;
phi=0.5;
chig=0.0039;
xi=0.1; %endo
lp=0.1; %endo
lg=0.1; %endo
%param production
ha=0.15;
hb=0.1;
alph=0.36;
sigmaa=0.1;%endo
sigmab=0.1;%endo
%param utility
zet=1.1;
gamm=0.8;

% ratio parameter
gi_Y=0.0724;
gc_Y=0.17;
b_Y=1.6;
T_Y=Ga*b_Y;
z_wp=0.26;

%consistency param
rhogc=0.9; rhogi=0.9; rhoz=0.9; rhotauc=0.9;
rhotaur=0.9; rhotauw=0.9; rhoMUg=0.9;
rhovg=0.9; rhoi=0.8; phiy=0.01;
phipi=1.1; rhoAf=0.9; rhoAm=0.9; rhoAg=0.9; rhochi=0.9;
rhonu=0.9;

varexo
epsgc epsgi epsz epstauw epstaur epstauc epsMUg epsvg epsi
epsAf epsAm epschi epsnu epsAg
;

model;
%1
np+ng+u=1;
%2
m=Am* v^kapp * u^(1-kapp);
%3
v=vp+vg;
%4
mu=m/v;
%5
fp=m/u*vp/v;
%6
fg=m/u*vg/v;
%7
u=(1-fg(-1)-fp(-1))*u(-1)+chi(-1)*np(-1)+chig*ng(-1);
%8
ng=(1-chig)*ng(-1)+fg(-1)*u(-1);
%9
Uc=(C/(1+gamm*nu*(zet-1)*n))^(-zet);
%10
Uc=(1+tauc)*lamb;
%11
lamb=psi*(1- d/2*(I/I(-1)-1)^2 - d*(I/I(-1)-1)*(I/I(-1)))
     +bet*psi(+1)*d*(I(+1)/I-1)*(I(+1)/I)^2;
%12
psi=bet*((1-taur)*lamb(+1)*r(+1)+psi(+1)*(1-delt));
%13
K=(1-delt)*K(-1)+(1-d/2*(I/I(-1)-1)^2)*I;
%14
N=z+bet*lamb(+1)/lamb*(fp*Op(+1)+fg*Og(+1)+(1-fp-fg)*N(+1));
%15
Op=(1-tauw)*wp*hp-nu*sigmaa*log(1-hp)/lamb
    +bet*lamb(+1)/lamb*((1-chi)*Op(+1)+chi*N(+1));
%16
Og=(1-tauw)*wg*hg-nu*sigmab*log(1-hg)/lamb
    +bet*lamb(+1)/lamb*((1-chig)*Og(+1)+chig*N(+1));
%17
J=(1-alph)*pm*y-wp*hp+bet*lamb(+1)/lamb*(1-chi)*J(+1);
%18
0=-lp+bet*lamb(+1)/lamb*mu*J(+1);
%19
y=Af * k(-1)^alph * hp^(1-alph) * GK(-1)^ha * gs^hb;
%20
Y=np*y;
%21
K=np*k;
%22
r=alph*pm*y/k(-1);
%23
wpstar*hp=phi*((1-alph)*pm*y+fp/mu*lp)
   +(1-phi)/(1-tauw)*(z-nu*sigmaa*log(1-hp)/lamb
   +bet*fg*lamb(+1)/lamb*(Og(+1)-N(+1)));
%24
pistar=et/(et-1)*pi*F1/F2;
%25
F1=lamb*pm*Y+bet*omegap*F1(+1)*pi(+1)^et;
%26
F2=lamb*Y+bet*omegap*F2(+1)*pi(+1)^(et-1);
%27
pi^(1-et)=omegap+(1-omegap)*pistar^(1-et);
%28
wp=omegaw*wp(-1)+(1-omegaw)*wpstar;
%29
b=i(-1)/pi*b(-1)+z*u+gc+gi+wg*ng*hg+lg*vg-T
  -tauw*wp*np*hp-tauw*wg*ng*hg-taur*r*K-tauc*C;
%30
log(gc)=rhogc*log(gc(-1))+(1-rhogc)*log(steady_state(gc))+epsgc;
%31
log(gi)=rhogi*log(gi(-1))+(1-rhogi)*log(steady_state(gi))+epsgi;
%32
log(z)=rhoz*log(z(-1))+(1-rhoz)*log(steady_state(z))+epsz;
%33
log(tauw)=rhotauw*log(tauw(-1))+(1-rhotauw)*log(steady_state(tauw))-epstauw;
%34
log(taur)=rhotaur*log(taur(-1))+(1-rhotaur)*log(steady_state(taur))-epstaur;
%35
log(tauc)=rhotauc*log(tauc(-1))+(1-rhotauc)*log(steady_state(tauc))-epstauc;
%36
gs=Ag*ng*hg;
%37
wg=(1+MUg)*wp;
%38
log(MUg)=rhoMUg*log(MUg(-1))+(1-rhoMUg)*log(steady_state(MUg))+epsMUg;
%39
log(vg)=rhovg*log(vg(-1))+(1-rhovg)*log(steady_state(vg))+epsvg;
%40
GK=(1-delt)*GK(-1)+gi;
%41
log(i)=rhoi*log(i(-1))+(1-rhoi)*log(steady_state(i))+phipi*log(pi)+phiy*log(y/steady_state(y))+epsi;
%42
lamb=bet*lamb(+1)*i/pi(+1);
%43
Y=C+I+gi+gc+u*z+vp*lp+vg*lg;
%44
T=Ga*b;
%45
log(Af)=rhoAf*log(Af(-1))+epsAf;
%46
log(Am)=rhoAm*log(Am(-1))+(1-rhoAm)*log(steady_state(Am))+epsAm;
%47
log(chi)=rhochi*log(chi(-1))+(1-rhochi)*log(steady_state(chi))+epschi;
%48
log(nu)=rhonu*log(nu(-1))+epsnu;
%49
n=np+ng;
%50
log(Ag)=rhoAg*log(Ag(-1))+(1-rhoAg)*log(steady_state(Ag))+epsAg;
%51
thet=v/u;
%52
Op=Og+xi;
%53
nu/(1-tauw)*sigmaa/lamb/(1-hp)=(1-alph)*pm*y/hp;
end;


steady;
check;

shocks;
var epstauc=0.01;
end;
stoch_simul(order=1,irf=100,periods=500);
