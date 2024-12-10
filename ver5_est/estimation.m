%% estimation_hpfilter data_construct
clc,clear;
cd 'C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_est'
data=readtable('.\data\data_after_clear.xls','VariableNamingRule','preserve');
%delete missing values (1985Q1-2019Q4)
data=data(23:160,:);

% seasonal adjustment
data.y=season_adj(data.y,'1985Q3');
data.c=season_adj(data.c,'1985Q3');
data.I=season_adj(data.I,'1985Q3');
data.tauc=season_adj(data.tauc,'1985Q3');
data.tauc=season_adj(data.tauc,'1985Q3');
data.tauk=season_adj(data.tauk,'1985Q3');
data.taun=season_adj(data.taun,'1985Q3');
data.z=season_adj(data.z,'1985Q3');
data.xi=season_adj(data.xi,'1985Q3');
data.ng=season_adj(data.ng,'1985Q3');
data.w=season_adj(data.w,'1985Q3');
data.b=season_adj(data.b,'1985Q3');
data.GC=season_adj(data.GC,'1985Q3');
data.GI=season_adj(data.GI,'1985Q3');
data.hp=season_adj(data.hour,'1985Q3');
data.v=season_adj(data.V,'1985Q3');

%one_side_hp_filter
data.v_obs=log(data.v)-one_sided_hp_filter_kalman(log(data.v),1600);
data.y_obs=log(data.y)-one_sided_hp_filter_kalman(log(data.y),1600);
data.c_obs=log(data.c)-one_sided_hp_filter_kalman(log(data.c),1600);
data.I_obs=log(data.I)-one_sided_hp_filter_kalman(log(data.I),1600);
data.w_obs=log(data.w)-one_sided_hp_filter_kalman(log(data.w),1600);
data.GC_obs=log(data.GC)-one_sided_hp_filter_kalman(log(data.GC),1600);
data.GI_obs=log(data.GI)-one_sided_hp_filter_kalman(log(data.GI),1600);
data.b_obs=log(data.b_pc)-one_sided_hp_filter_kalman(log(data.b),1600);
data.z_obs=log(data.z)-one_sided_hp_filter_kalman(log(data.z),1600);
data.tauc_obs=log(data.tauc)-one_sided_hp_filter_kalman(log(data.tauc),1600);
data.tauk_obs=log(data.tauk)-one_sided_hp_filter_kalman(log(data.tauk),1600);
data.taun_obs=log(data.taun)-one_sided_hp_filter_kalman(log(data.taun),1600);
data.xi_obs=log(data.xi)-one_sided_hp_filter_kalman(log(data.xi),1600);
data.i_obs=log(data.i)-one_sided_hp_filter_kalman(log(data.i),1600);
data.pi_obs=log(data.pi)-one_sided_hp_filter_kalman(log(data.pi),1600);
data.hp_obs=log(data.hp)-one_sided_hp_filter_kalman(log(data.hp),1600);
data.theta_obs=log(data.theta)-one_sided_hp_filter_kalman(log(data.theta),1600);

data.n_ob=log(data.n)-one_sided_hp_filter_kalman(log(data.n),1600);
data.ng_obs=log(data.ng)-one_sided_hp_filter_kalman(log(data.ng),1600);
data.u_obs=log(data.u)-one_sided_hp_filter_kalman(log(data.u),1600);


data.v=season_adj(data.n.*data.theta,'1985Q3');
data.v_obs=log(data.v)-one_sided_hp_filter_kalman(log(data.v),1600);

%obs_variables
uobs=data.u_obs;
vobs=data.v_obs;
thetaobs=data.theta_obs;
yobs=data.y_obs;
cobs=data.c_obs;
Iobs=data.I_obs;
wobs=data.w_obs;
piobs=data.pi_obs;
gcobs=data.GC_obs;
giobs=data.GI_obs;
zobs=data.z_obs;
taucobs=data.tauc_obs;
tauwobs=data.taun_obs;
taurobs=data.tauk_obs;
ngobs=data.ng_obs;
iobs=data.i_obs;
nob=data.n_ob;
hpobs=data.hp_obs;

cd 'C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_est'
save('dataset_hp.mat','yobs','cobs','Iobs','wobs','piobs',...
    'gcobs','giobs','zobs','ngobs','iobs','nob',...
    'taucobs','tauwobs','taurobs','hpobs','thetaobs',...
    'uobs','vobs');





%% estimation start
clc,clear;
dynare est_search_gov_emp
est_oo_=oo_;
est_M_=M_;
save('est_result','est_oo_','est_M_');