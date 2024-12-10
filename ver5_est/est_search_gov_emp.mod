var
%main
np ng u m Am v vp vg mu fp
fg chi nu n tauc lamb psis I r Y
K N z Op Og wp tauw wg J pm
y GK wpstar pistar pi F1 F2 b gc k
gi i T C taur gs Ag MUg Af Uc 
thet hp hg et
% observation datas
yobs cobs Iobs wobs piobs gcobs giobs
bobs zobs ngobs iobs nob taucobs taurobs tauwobs hpobs thetaobs
% additional
creve wreve rreve %tax revenues
vgeve             %government employment revenues
benefit %unemployment benefit
wupper wlower %reservation wages
grevenue %government revenue
gexpend %government expenditure
deltb   %change in government debt
;

parameters
kapp chig gamm zet d bet delt alph lp lg ets
omegaw omegap sigmaa sigmab xi ha hb
rhogc rhogi rhoz rhotauw rhotauc rhotaur 
phi rhoMUg rhovg rhoi phipi phiy Ga rhoAf rhoAm rhoAg
rhonu rhochi rhoet gi_Y gc_Y b_Y T_Y z_wp
taucs taurs tauws MUgs nus ngs hgs hps us thets
gss chis

phibtauc phiytauc
phibtaur phiytaur
phibtauw phiytauw
phibgc phiygc
phibgi phiygi
;
% basic param
delt=0.025;
bet=0.96;
d=0.2472;
% price param
omegap=0.6190;
omegaw=0.5956;
ets=11;
%param wage and match
kapp=0.5025;
phi=0.6691;
chis=0.012;
chig=0.012;
nus=1;
xi=0.1; %endo 
lp=0.1; %endo
lg=0.1; %endo
%param production
ha=0.1567; hb=0.1016;
alph=0.36;
sigmaa=0.1;%endo
sigmab=0.1;%endo

%param utility
zet=1.2147; gamm=0.3388;

% ratio parameter
ngs=0.07; gss=1/3;

hgs=0.381; us=0.034;
hps=0.381; thets=1.0;

Ga=0.06; % lump-sum tax/debt
gi_Y=0.0724; gc_Y=0.17; b_Y=1.6; 
T_Y=Ga*b_Y; z_wp=0.26;
MUgs=0.127; % markup rate of gov wage

%taxes
tauws=0.26; taucs=0.07; taurs=0.44;

%consistency param
rhoAm=0.6801; 
rhochi=0.8605;
rhonu=0.7526; 
rhoAf=0.9882; 
rhoet=0.8;
rhogc=0.7710; 
rhogi=0.8310; 
rhovg=0.7646; 
rhoMUg=0.7843;
rhotauc=0.9883;
rhotaur=0.7533; 
rhotauw=0.8; 
rhoi=0.7581; 
rhoz=0.7589; 
rhoAg=0.7369;

phipi=1.0465; 
phiy=0.0441;
%fiscal rules param
phibgc=-0.0904; phiygc=-0.1026;
phibgi=-0.1006; phiygi=-0.1171;
phibtauc=0.0783;  phiytauc=0.0532;
phibtaur=0.0473; phiytaur=0.0420;
phibtauw=0.0096; phiytauw=0.0646;


varexo
epsgc epsgi epsz epstauw epstaur epstauc 
epsMUg epsvg
%epsAg
epsi epsAf epsAm epschi epsnu epset
% observation
ey ez ew
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
Uc=(C/(1+nu*gamm*(zet-1)*n))^(-zet);
%10
Uc=(1+tauc)*lamb;
%11
lamb=psis*(1- d/2*(I/I(-1)-1)^2 - d*(I/I(-1)-1)*(I/I(-1)))
     +bet*psis(+1)*d*(I(+1)/I-1)*(I(+1)/I)^2;
%12
psis=bet*((1-taur)*lamb(+1)*r(+1)+psis(+1)*(1-delt));
%13
K=(1-delt)*K(-1)+(1-d/2*(I/I(-1)-1)^2)*I;
%14
N=z+bet*lamb(+1)/lamb*(fp*Op(+1)+fg*Og(+1)+(1-fp-fg)*N(+1));
%15
Op=(1-tauw)*wp*hp+sigmaa*log(1-hp)/lamb
    +bet*lamb(+1)/lamb*((1-chi)*Op(+1)+chi*N(+1));
%16
Og=(1-tauw)*wg*hg+sigmab*log(1-hg)/lamb
    +bet*lamb(+1)/lamb*((1-chig)*Og(+1)+chig*N(+1));
%17
%J=(1-alph)*pm*y-wp*hp+bet*lamb(+1)/lamb*(1-chi)*J(+1);
J=(1-alph)*pm*y-wp*hp;
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
% wpstar*hp=phi*((1-alph)*pm*y+fp/mu*lp)
%    +(1-phi)/(1-tauw)*(z-nu*sigmaa*log(1-hp)/lamb
%    +bet*fg*lamb(+1)/lamb*(Og(+1)-N(+1)));
wpstar*hp=phi*((1-alph)*pm*y+fp/mu*lp)
   +(1-phi)/(1-tauw)*(z-sigmaa*log(1-hp)/lamb
   +bet*fg*lamb(+1)/lamb*(Og(+1)-N(+1)))
   +bet*(1-chi-fp)*lamb(+1)/lamb*(wpstar*hp-wp*hp);
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
log(gc)=rhogc*log(gc(-1))+(1-rhogc)*(log(steady_state(gc))
       +phibgc*log(b/steady_state(b))
       +phiygc*log(y/steady_state(y)))+epsgc;
%31
log(gi)=rhogi*log(gi(-1))+(1-rhogi)*(log(steady_state(gi))
        +phibgi*log(b/steady_state(b))
        +phiygi*log(y/steady_state(y)))+epsgi;
%32
log(z)=rhoz*log(z(-1))+(1-rhoz)*log(steady_state(z))+epsz;
%33
log(tauw)=rhotauw*log(tauw(-1))+(1-rhotauw)*(log(steady_state(tauw))
        +phibtauw*log(b/steady_state(b))
        +phiytauw*log(y/steady_state(y)))+epstauw;
%34
log(taur)=rhotaur*log(taur(-1))+(1-rhotaur)*(log(steady_state(taur))
        +phibtaur*log(b/steady_state(b))
        +phiytaur*log(y/steady_state(y)))+epstaur;
%35
log(tauc)=rhotauc*log(tauc(-1))+(1-rhotauc)*(log(steady_state(tauc))
        +phibtauc*log(b/steady_state(b))
        +phiytauc*log(y/steady_state(y)))+epstauc;
%36
gs=Ag*ng*hg;
%37
wg=(1+MUg)*wp;
%38
log(MUg)=rhoMUg*log(MUg(-1))+(1-rhoMUg)*log(steady_state(MUg))+epsMUg;
%log(MUg)=rhoMUg*log(MUg(-1))+(1-rhoMUg)*log(steady_state(MUg));
%39
log(vg)=rhovg*log(vg(-1))+(1-rhovg)*log(steady_state(vg))+epsvg;
%log(vg)=rhovg*log(vg(-1))+(1-rhovg)*log(steady_state(vg));

%40
GK=(1-delt)*GK(-1)+gi;
%41
log(i)=rhoi*log(i(-1))+(1-rhoi)*(log(steady_state(i))+phipi*log(pi)+phiy*log(y/steady_state(y)))+epsi;
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
log(nu)=rhonu*log(nu(-1))+(1-rhonu)*log(steady_state(nu))+epsnu;
%49
n=np+ng;
%50
%log(Ag)=rhoAg*log(Ag(-1))+(1-rhoAg)*log(steady_state(Ag))+epsAg;
log(Ag)=rhoAg*log(Ag(-1))+(1-rhoAg)*log(steady_state(Ag));
%51
thet=v/u;
%52
Op=Og+xi;
%53
nu/(1-tauw)*sigmaa/lamb/(1-hp)=(1-alph)*pm*y/hp;
%54
log(et)=rhoet*log(et(-1))+(1-rhoet)*log(steady_state(et))+epset;
%------------------------------
%observation equations
yobs=log(Y)-log(steady_state(Y))+ey;
cobs=log(C)-log(steady_state(C));
Iobs=log(I)-log(steady_state(I));
wobs=log(wp)-log(steady_state(wp))+ew;
piobs=log(pi)-log(steady_state(pi));
gcobs=log(gc)-log(steady_state(gc));
giobs=log(gi)-log(steady_state(gi));
bobs=log(b)-log(steady_state(b));
zobs=log(z)-log(steady_state(z))+ez;
taucobs=log(tauc)-log(steady_state(tauc));
tauwobs=log(tauw)-log(steady_state(tauw));
taurobs=log(taur)-log(steady_state(taur));
iobs=log(i)-log(steady_state(i));
hpobs=log(hp)-log(steady_state(hp));
thetaobs=log(thet)-log(steady_state(thet));

ngobs=log(ng)-log(steady_state(ng));
nob=log(n)-log(steady_state(n));

%--------------------------------
% additional
%tax revenues
creve=tauc*C;
rreve=taur*r*K;
wreve=tauw*wp*np*hp+tauw*wg*ng*hg;
vgeve=lg*vg+wg*ng*hg;
% unemployment benefit
benefit=z*u;
%researvation wages
wupper=(1-alph)*pm*y+fp/mu*lp;
wlower=z-nu*sigmaa*log(1-hp)/lamb
   +bet*fg*lamb(+1)/lamb*(Og(+1)-N(+1));
%government revenue and expenditure
grevenue=T+tauw*wp*np*hp+tauw*wg*ng*hg+taur*r*K+tauc*C;
gexpend=z*u+gc+gi+wg*ng*hg+lg*vg;
deltb=b-b(-1);
end;

steady;
check;

shocks;
var epstauc=0.01;
var epstaur=0.01;
var epstauw=0.01;
var epsgc=0.01;
var epsgi=0.01;
var epsi=0.01;
var epsz=0.01;
var epsAf=0.01;
var epsnu=0.01;
var epschi=0.01;
var epsAm=0.01;
var epset=0.01;
%--------------------
var ey=0.01^2;
var ew=0.01^2;
var ez=0.01^2;
end;

varobs
yobs cobs wobs piobs gcobs giobs iobs
zobs ngobs nob thetaobs
taucobs taurobs tauwobs hpobs
;

%  estimated_params_init(use_calibration);
%  end;

 estimated_params;
%%param production
%ha,beta_pdf,0.1,0.03;
%hb,beta_pdf,0.1,0.03;

%param wage and match
%kapp,beta_pdf,0.6,0.15;
%phi,beta_pdf,0.5,0.15;

%%param utility
% zet,gamma_pdf,1.2,0.5; 
% gamm,beta_pdf,0.4,0.15;

%%basic param
% price param
 omegap,beta_pdf,0.6,0.15;
%  omegaw,beta_pdf,0.6,0.15;
 d,beta_pdf,0.3,0.15;


% %consistency param
 rhoAm,beta_pdf,0.8,0.15;
 rhochi,beta_pdf,0.8,0.15;
 rhonu,beta_pdf,0.8,0.15;
 rhoAf,beta_pdf,0.8,0.15;
 rhoet,beta_pdf,0.8,0.15;
 rhogc,beta_pdf,0.8,0.15;
 rhogi,beta_pdf,0.8,0.15;
 rhovg,beta_pdf,0.8,0.15;
 rhoMUg,beta_pdf,0.8,0.15;
 rhotauc,beta_pdf,0.8,0.15;
 rhotaur,beta_pdf,0.8,0.15;
 rhotauw,beta_pdf,0.8,0.15;
 rhoi,beta_pdf,0.8,0.15;
 rhoz,beta_pdf,0.8,0.15;

% fiscal rule param
% phibtauc,normal_pdf,0,0.05;
% phibtauw,normal_pdf,0,0.05;
% phibtaur,normal_pdf,0,0.05;
% phibgc,normal_pdf,0,0.05;
% phiygi,normal_pdf,0,0.05;

% phiytauc,normal_pdf,0,0.05;
% phiytauw,normal_pdf,0,0.05;
% phiytaur,normal_pdf,0,0.05;
% phiygc,normal_pdf,0,0.05;
% phibgi,normal_pdf,0,0.05;

%interest rate parameter
% phipi,gamma_pdf,1.1,0.05;
% phiy,gamma_pdf,0.1,0.05;
 
% %shocks (14)
stderr epsAf,inv_gamma_pdf,0.5,Inf;
stderr epsAm,inv_gamma_pdf,0.5,Inf;
stderr epsnu,inv_gamma_pdf,0.5,Inf;
stderr epset,inv_gamma_pdf,0.5,Inf;
stderr epschi,inv_gamma_pdf,0.5,Inf;
stderr epsgc,inv_gamma_pdf,0.5,Inf;
stderr epsgi,inv_gamma_pdf,0.5,Inf;
stderr epsvg,inv_gamma_pdf,0.5,Inf;
stderr epsMUg,inv_gamma_pdf,0.5,Inf;
stderr epstauc,inv_gamma_pdf,0.5,Inf;
stderr epstaur,inv_gamma_pdf,0.5,Inf;
stderr epstauw,inv_gamma_pdf,0.5,Inf;
stderr epsi,inv_gamma_pdf,0.5,Inf;
stderr epsz,inv_gamma_pdf,0.5,Inf;

% observation shock
stderr ey,inv_gamma_pdf,0.15,0.3; 
stderr ew,inv_gamma_pdf,0.15,0.3;
stderr ez,inv_gamma_pdf,0.15,0.3;
end;

estimation(order=1,datafile='dataset_hp.mat',mode_check,mh_replic=400000,mode_compute=6,mh_nblocks=2,mh_drop=0.4,mh_jscale=0.8,prior_trunc=0,nograph);
shock_decomposition(datafile='dataset_hp.mat',nograph) yobs nob piobs;
stoch_simul(order=1,irf=20,periods=0,nograph);
plot_shock_decomposition(write_xls) yobs nob piobs;
