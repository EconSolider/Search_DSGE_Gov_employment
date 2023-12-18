function [ys,params,check] = e_search_gov_emp_steadystate(ys,exo,M_,options_)
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] 0 if steady state computation worked and to
%                        1 of not (allows to impose restrictions on parameters)

%% Step 0: initialize indicator and set options for numerical solver
check = 0;
options = optimset('Display','Final','TolX',1e-10,'TolFun',1e-10);

%% Step: read out parameters to access them with their name
for ii = 1:M_.param_nbr
  eval([ M_.param_names{ii} ' = M_.params(' int2str(ii) ');']);
end
%% 
ng=0.1;
MUg=0.127;
tauw=0.26;
tauc=0.07;
taur=0.44;
chi=0.012;
pi=1;
pistar=1;
pm=(et-1)/et;
r=(1/bet+delt-1)/(1-taur);
K_Y=alph*pm/r;
k_y=K_Y;
I_Y=delt*K_Y;
GK_Y=1/delt*gi_Y;
%labor market
u=0.026;
hp=1/3;
thet=1.1;
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
gs=1/3;
y=(k_y^alph*hp^(1-alph)*(GK_Y*np)^ha*gs^hb)^(1/(1-alph-ha));
Y=np*y;
I=I_Y*Y;
K=K_Y*Y;
k=K/np;
gi=gi_Y*Y;
GK=GK_Y*Y;
b=b_Y*Y;
T=T_Y*Y;
hg=1/3;
Ag=gs/ng/hg;
gc=gc_Y*Y;

% c w lp lg
fun_a=@(x)Fun_a(x,alph,pm,y,hp,chi,mu,bet,I,gi,gc,u,z_wp,vp,vg,Y,b,i,T,ng,tauw,MUg,np,hg,taur,r,K,tauc);
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
sigmaa=(1-alph)*pm*y/hp*(1-tauw)*lamb*(1-hp);
% value functions
fun_b=@(x)Fun_b(x,z,bet,fp,fg,tauw,wp,hp,sigmaa,lamb,chi,wg,hg,chig,phi,alph,pm,y,mu,lp);
x0=[0.1 10 10 10 0.05];
x = fsolve(fun_b, x0, options);
N=x(1);
Op=x(2);
Og=x(3);
xi=x(4);
sigmab=x(5);

F1=lamb*pm*Y/(1-bet*omegap);
F2=lamb*Y/(1-bet*omegap);
nu=1;
psi=lamb;
wpstar=wp;
Af=1;


%% Step: Update parameters and variables
params=NaN(M_.param_nbr,1);
for iter = 1:M_.param_nbr %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

for ii = 1:M_.orig_endo_nbr %auxiliary variables are set automatically
  eval(['ys(' int2str(ii) ') = ' M_.endo_names{ii} ';']);
end

end