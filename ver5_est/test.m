%% parameters
clc,clear;
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
zet=0.2147; gamm=0.3388;

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
%% endovar
%%
et=ets;
ng=ngs;
MUg=MUgs;
tauw=tauws;
tauc=taucs;
taur=taurs;
chi=chis;
nu=nus;
pi=1;
pistar=1;
pm=(et-1)/et;
r=(1/bet+delt-1)/(1-taur);
K_Y=alph*pm/r;
k_y=K_Y;
I_Y=delt*K_Y;
GK_Y=1/delt*gi_Y;
%labor market
u=us;
hp=hps;
thet=thets;
v=thet*u;
np=1-ng-u;
fg=chig*ng/u;
fp=1-(u-chi*np-chig*ng)/u-fg;
vg=(1+fp/fg)^(-1)*v;
vp=v-vg;
m=fg*v*u/vg;
Am=m/(v^kapp*u^(1-kapp));
n=np+ng;
mu=m/v;
i=1/bet;
%production
gs=gss;
y=(k_y^alph*hp^(1-alph)*(GK_Y*np)^ha*gs^hb)^(1/(1-alph-ha));
Y=np*y;
I=I_Y*Y;
K=K_Y*Y;
k=K/np;
gi=gi_Y*Y;
GK=GK_Y*Y;
b=b_Y*Y;
T=T_Y*Y;
hg=hgs;
Ag=gs/ng/hg;
gc=gc_Y*Y;

% c w lp lg
fun_a=@(x)Fun_a(x,alph,pm,y,hp,mu,bet,I,gi,gc,u,z_wp,vp,vg,Y,b,i,T,ng,tauw,MUg,np,hg,taur,r,K,tauc);
x0=[0.1 0.1 0.1];
options = optimoptions('fsolve', 'Display', 'none','TolX', 1e-8,'TolFun', 1e-8);
x = fsolve(fun_a, x0, options);
C=x(1);
lp=x(2);
wp=x(3);
lg=0.15*lp;

if C<0
    error('C<0')
end

if lp<0
    error('lp<0')
end

if wp<0
    error('wp<0')
end

z=z_wp*wp;
J=lp/bet/mu;
Uc=(C/(1+(zet-1)*gamm*n))^(-zet);
lamb=(C/(1+(zet-1)*gamm*n))^(-zet)/(1+tauc);
wg=(1+MUg)*wp;
sigmaa=(1-alph)*pm*y/hp*(1-tauw)*lamb*(1-hp)/nu;
% value functions
fun_b=@(x)Fun_b(x,z,bet,fp,fg,tauw,wp,hp,sigmaa,lamb,chi,wg,hg,chig,phi,alph,pm,y,mu,lp,nu);
x0=[0.1 10 10 10 0.05];
x = fsolve(fun_b, x0, options);
N=x(1);
Op=x(2);
Og=x(3);
xi=x(4);
sigmab=x(5);

% if sigmab<0
% error('sigmab<0')
% end


F1=lamb*pm*Y/(1-bet*omegap);
F2=lamb*Y/(1-bet*omegap);

psis=lamb;
wpstar=wp;
Af=1;