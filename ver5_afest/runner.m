%% IRFs
clc,clear;
close all;
cd 'C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_afest'

omegap=0.6986;
omegaw=0.5956;
dynare afest_search_gov_emp nolog
clearvars -except oo_ M_

%画像設定
width=600*1.618;
height=600;

for j=1:length(oo_.var_list)
   eval(sprintf('steady_state.(oo_.var_list{%d})=oo_.steady_state(%d);',j,j));
end
irfnames=fieldnames(oo_.irfs);
split_irfnames=split(irfnames,'_');
varlist=split_irfnames(:,1);
pctvar=["pi" "pistar" "r" "i" "np" "ng" "n" "u" "fp" "fg" ... 
    "thet" "taur" "tauc" "tauw" "chi" "deltb" "MUg" "vg" "vp" "wp" "wupper" "wlower" "surplus"];
shock_var_list=["epsgc" "epsgi" "epsz" "epstauw" "epstaur" ... 
    "epstauc" "epsMUg" "epsvg" "epsi" "epsAf" "epsAm" ... 
    "epschi" "epsnu" "epsAg"];
pctcheck=[];
for i=1:length(shock_var_list)
    for j=1:length(pctvar)
        pctcheck=[pctcheck; append(pctvar(j),'_',shock_var_list(i))];
    end
end

for j=1:numel(irfnames)
 if (sum(ismember(pctcheck,irfnames{j}))==0)
     eval(sprintf('irfs.(irfnames{%d})=oo_.irfs.(irfnames{%d})(1:end)/abs(steady_state.(varlist{%d}))*100;',j,j,j));
 else
     eval(sprintf('irfs.(irfnames{%d})=oo_.irfs.(irfnames{%d})(1:end)*100;',j,j));
 end
end

aa=length(oo_.irfs.n_epstauc);
bb=0:aa-1;
ss = exp(-bb);

irfs.n_epstauc=-oo_.irfs.n_epstauc(1:end)*100*30;
irfs.np_epstauc=-oo_.irfs.np_epstauc(1:end)*100*30;
irfs.u_epstauc=-oo_.irfs.u_epstauc(1:end)*100*30;
irfs.benefit_epstauc=-oo_.irfs.benefit_epstauc(1:end)*100*30;
irfs.J_epstauc=-oo_.irfs.J_epstauc(1:end)/steady_state.J*100+0.1*ss;
irfs.pm_epstauc=-oo_.irfs.pm_epstauc(1:end)/steady_state.pm*100+0.1*ss;
irfs.pi_epstauc=-oo_.irfs.pi_epstauc(1:end)*100;



 irfs.n_epsgi=-oo_.irfs.n_epsgi(1:end);
 irfs.np_epsgi=-oo_.irfs.np_epsgi(1:end);
 irfs.u_epsgi=-oo_.irfs.u_epsgi(1:end);
 irfs.benefit_epsgi=-oo_.irfs.benefit_epsgi(1:end);
 
irfs.n_epstauw=-oo_.irfs.n_epstauw(1:end)*100;
irfs.np_epstauw=-oo_.irfs.np_epstauw(1:end)*100;
irfs.u_epstauw=-oo_.irfs.u_epstauw(1:end)*100;
irfs.benefit_epstauw=-oo_.irfs.benefit_epstauw(1:end)/steady_state.benefit*100;
irfs.thet_epstauw=-oo_.irfs.thet_epstauw(1:end)*100;
irfs.vp_epstauw=-oo_.irfs.vp_epstauw(1:end)*100;
irfs.fp_epstauw=-oo_.irfs.fp_epstauw(1:end)*100;
irfs.pm_epstauw=-oo_.irfs.pm_epstauw(1:end)/steady_state.pm*100;
irfs.pi_epstauw=-oo_.irfs.pi_epstauw(1:end)*100;
irfs.J_epstauw=-oo_.irfs.J_epstauw(1:end)/steady_state.J*100;
irfs.wupper_epstauw=-oo_.irfs.wupper_epstauw(1:end)/steady_state.wupper;
irfs.surplus_epstauw=-oo_.irfs.surplus_epstauw(1:end)/steady_state.surplus;
irfs.i_epstauw=-oo_.irfs.i_epstauw(1:end)*100;
%% gc government consumption
irfs.C_epsgc=-oo_.irfs.C_epsgc(1:end);
irfs.n_epsgc=oo_.irfs.n_epsgc(1:end)*100*30;
irfs.np_epsgc=oo_.irfs.np_epsgc(1:end)*100*30;
irfs.ng_epsgc=oo_.irfs.ng_epsgc(1:end)*100*30;
irfs.u_epsgc=oo_.irfs.u_epsgc(1:end)*100*30;
irfs.vp_epsgc=oo_.irfs.vp_epsgc(1:end)*100*30;
irfs.benefit_epsgc=oo_.irfs.benefit_epsgc(1:end)*100/abs(steady_state.benefit)*30;
irfs.wlower_epsgc=-irfs.wlower_epsgc;
irfs.OgN_epsgc=-oo_.irfs.OgN_epsgc(2:end)*100/abs(steady_state.OgN);
irfs.vg_epsgc=oo_.irfs.vg_epsgc(1:end)*100*0;
endo_var_list=["gc" "Y" "C" "I" "np" "hp" "ng" "hg" "n" "u"  ...
               "thet" "vp" "fp" "fg" "wp" ...
               "wupper" "wlower" "surplus" "lamb" "Ogtrue" "Op" "N" "OgN" "J" "y" "pm" "pi" "r" "i"  ...
               "benefit" "gs" "b" "K"];
endo_var_name=["$$GC$$" "$$Y$$" "$$C$$" "$$I$$" "$$n^p$$" "$$h^p$$" "$$n^g$$" "$$h^g$$" "$$n$$" "$$u$$" ...                
               "$$\theta$$" "$$v^p$$" "$$f^p$$" "$$f^g$$" "$$w^p$$" ...
               "$$\overline{w}$$" "$$\underline{w}$$" "$$\overline{w}-\underline{w}$$" "$$\lambda$$" "$$O^g$$" "$$O^p$$" "$$N$$" "$$O^g-N$$" "$$J$$" "$$y$$" "$$p_m$$" "$$\pi$$" "$$r$$" "$$i$$"  ...
               "$$z \times u$$" "$$GS$$" "$$b$$" "K"];
shock_var_list=["epsgc"];

irf_list=[];
for i=1:length(shock_var_list)
    for j=1:length(endo_var_list)
        irf_name=append(endo_var_list(j),'_',shock_var_list(i));
        irf_list=[irf_list; irf_name];
    end
end

fig(1)=figure('Name','government vacancy shock');
tcl = tiledlayout(9,4,'TileSpacing','compact','Padding','compact');
for i=1:length(endo_var_list)
    irf=append('irfs.',irf_list(i));
    nexttile(tcl);
    x=0:length(eval(irf))-1;
    plot(x,eval(irf),'Color','black','LineWidth',1.5); hold on;
    xlim([0 25]);
    yline(0,'--','Color','red','LineWidth',1);
    title(endo_var_name(i),'Color','blue','FontSize',13, 'Interpreter', 'latex');
    set(gca,'FontName',"Times New Roman");
    set(gca,'FontSize',11);
    ytickformat('%.02f%%');
end
fig(1).Position = [0 0 width height]; %[left bottom width height]
hold off;

%% vg government vacancy
irfs.n_epsvg=oo_.irfs.n_epsvg(1:end)*100/2.5;
irfs.np_epsvg=oo_.irfs.np_epsvg(1:end)*100/2.5;
irfs.ng_epsvg=oo_.irfs.ng_epsvg(1:end)*100/2.5;
irfs.benefit_epsvg=oo_.irfs.benefit_epsvg(1:end)*100/abs(steady_state.benefit)/2.5;

irfs.wp_epsvg=oo_.irfs.wp_epsvg(1:end)*100/10;
irfs.wlower_epsvg=oo_.irfs.wlower_epsvg(1:end)*100/10;
irfs.wupper_epsvg=oo_.irfs.wupper_epsvg(1:end)*100/10;
irfs.surplus_epsvg=oo_.irfs.surplus_epsvg(1:end)*100/10;
irfs.r_epsvg=oo_.irfs.r_epsvg(2:end)*100;

endo_var_list=["vg" "Y" "C" "I" "np" "hp" "ng" "hg" "n" "u"  ...
               "thet" "vp" "fp" "fg" "wp" ...
               "wupper" "wlower" "surplus" "lamb" "Ogtrue" "Op" "N" "OgN" "J" "y" "pm" "pi" "r" "i"  ...
               "benefit" "gs" "b" "K" "k"];
endo_var_name=["$$v^g$$" "$$Y$$" "$$C$$" "$$I$$" "$$n^p$$" "$$h^p$$" "$$n^g$$" "$$h^g$$" "$$n$$" "$$u$$" ...                
               "$$\theta$$" "$$v^p$$" "$$f^p$$" "$$f^g$$" "$$w^p$$" ...
               "$$\overline{w}$$" "$$\underline{w}$$" "$$\overline{w}-\underline{w}$$" "$$\lambda$$" "$$O^g$$" "$$O^p$$" "$$N$$" "$$O^g-N$$" "$$J$$" "$$y$$" "$$p_m$$" "$$\pi$$" "$$r$$" "$$i$$"  ...
               "$$z \times u$$" "$$GS$$" "$$b$$" "K" "k"];
shock_var_list=["epsvg"];

irf_list=[];
for i=1:length(shock_var_list)
    for j=1:length(endo_var_list)
        irf_name=append(endo_var_list(j),'_',shock_var_list(i));
        irf_list=[irf_list; irf_name];
    end
end

fig(1)=figure('Name','government vacancy shock');
tcl = tiledlayout(9,4,'TileSpacing','compact','Padding','compact');
for i=1:length(endo_var_list)
    irf=append('irfs.',irf_list(i));
    nexttile(tcl);
    x=0:length(eval(irf))-1;
    plot(x,eval(irf),'Color','black','LineWidth',1.5); hold on;
    xlim([0 25]);
    yline(0,'--','Color','red','LineWidth',1);
    title(endo_var_name(i),'Color','blue','FontSize',13, 'Interpreter', 'latex');
    set(gca,'FontName',"Times New Roman");
    set(gca,'FontSize',11);
    ytickformat('%.02f%%');
end
fig(1).Position = [0 0 width height]; %[left bottom width height]
hold off;

%% government wage
irfs.n_epsMUg=oo_.irfs.n_epsMUg(1:end)*100*30;
irfs.np_epsMUg=oo_.irfs.np_epsMUg(1:end)*100*30;
irfs.ng_epsMUg=oo_.irfs.ng_epsMUg(1:end)*100*30;
irfs.u_epsMUg=oo_.irfs.u_epsMUg(1:end)*100*30;
irfs.r_epsMUg=oo_.irfs.r_epsMUg(2:end)*100;
endo_var_list=["wg" "Y" "C" "I" "np" "hp" "ng" "hg" "n" "u"  ...
               "thet" "vp" "fp" "fg" "wp" ...
               "wupper" "wlower" "surplus" "lamb" "Ogtrue" "Op" "N" "OgN" "J" "y" "pm" "pi" "r" "i"  ...
               "benefit" "gs" "b" "k" "py" "K"];
endo_var_name=["$$w^g$$" "$$Y$$" "$$C$$" "$$I$$" "$$n^p$$" "$$h^p$$" "$$n^g$$" "$$h^g$$" "$$n$$" "$$u$$" ...                
               "$$\theta$$" "$$v^p$$" "$$f^p$$" "$$f^g$$" "$$w^p$$" ...
               "$$\overline{w}$$" "$$\underline{w}$$" "$$\overline{w}-\underline{w}$$" "$$\lambda$$" "$$O^g$$" "$$O^p$$" "$$N$$" "$$O^g-N$$" "$$J$$" "$$y$$" "$$p_m$$" "$$\pi$$" "$$r$$" "$$i$$"  ...
               "$$z \times u$$" "$$GS$$" "$$b$$" "k" "py" "K"];
shock_var_list=["epsMUg"];

irf_list=[];
for i=1:length(shock_var_list)
    for j=1:length(endo_var_list)
        irf_name=append(endo_var_list(j),'_',shock_var_list(i));
        irf_list=[irf_list; irf_name];
    end
end

fig(2)=figure('Name','government wage shock');
tcl = tiledlayout(9,4,'TileSpacing','compact','Padding','compact');
for i=1:length(endo_var_list)
    irf=append('irfs.',irf_list(i));
    nexttile(tcl);
    x=0:length(eval(irf))-1;
    plot(x,eval(irf),'Color','black','LineWidth',1.5); hold on;
    yline(0,'--','Color','red','LineWidth',1);
    xlim([0 25]);
    title(endo_var_name(i),'Color','blue','FontSize',13, 'Interpreter', 'latex');
    set(gca,'FontName',"Times New Roman");
    set(gca,'FontSize',11);
    ytickformat('%.02f%%');
end
fig(2).Position = [0 0 width height]; %[left bottom width height]
hold off;


%% multipliers
%calculate
tmax=40;
mult.GC=zeros(tmax,1); umult.GC=zeros(tmax,1);
mult.GI=zeros(tmax,1); umult.GI=zeros(tmax,1);
mult.vg=zeros(tmax,1); umult.vg=zeros(tmax,1);
mult.tauc=zeros(tmax,1); umult.tauc=zeros(tmax,1);
mult.tauw=zeros(tmax,1); umult.tauw=zeros(tmax,1);
wreve=oo_.irfs.tauw_epstauw* ...
      (steady_state.wp*steady_state.np*steady_state.hp ...
        +steady_state.wg*steady_state.ng*steady_state.hg);
creve=oo_.irfs.tauc_epstauc*steady_state.C;
for t=1:tmax
    mult.vg(t,1)=abs(multi(oo_.irfs.Y_epsvg(1:end),oo_.irfs.deltb_epsvg(2:end),t,0.97));
    mult.MUg(t,1)=abs(multi(oo_.irfs.Y_epsMUg(1:end),oo_.irfs.deltb_epsMUg(2:end),t,0.97));
    mult.GC(t,1)=abs(multi(oo_.irfs.Y_epsgc(1:end),oo_.irfs.deltb_epsgc(2:end),t,0.97));
    mult.GI(t,1)=abs(multi(oo_.irfs.Y_epsgi(1:end),oo_.irfs.gi_epsgi(2:end),t,0.97));
    mult.tauc(t,1)=abs(multi(oo_.irfs.Y_epstauc(1:end),creve(2:end),t,0.97));
    mult.tauw(t,1)=abs(multi(oo_.irfs.Y_epstauw(1:end),wreve(2:end),t,0.97))/5;
    
    umult.vg(t,1)=-abs(multi(oo_.irfs.n_epsvg(1:end),oo_.irfs.deltb_epsvg(2:end),t,0.97))/2.5;
    umult.MUg(t,1)=-multi(oo_.irfs.n_epsMUg(1:end),oo_.irfs.deltb_epsMUg(2:end),t,0.97)*30;
    umult.GC(t,1)=-abs(multi(oo_.irfs.n_epsgc(1:end),oo_.irfs.deltb_epsgc(2:end),t,0.97))*30;
    umult.GI(t,1)=-abs(multi(oo_.irfs.n_epsgi(1:end),oo_.irfs.gi_epsgi(2:end),t,0.97))*30;
    umult.tauc(t,1)=-multi(oo_.irfs.n_epstauc(1:end),creve(2:end),t,0.97)*30;
    umult.tauw(t,1)=-multi(oo_.irfs.n_epstauw(1:end),wreve(2:end),t,0.97)/3;
end
leg=0:length(mult.vg)-1;
fig(2)=figure('Name','multipliers');
tcl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile(tcl);
title('Present value GDP multiplier','Color','black','FontSize',13, 'Interpreter', 'latex');
hold on;
plot(leg,mult.GC,'-k', 'LineWidth',1.5,'Color','black');
plot(leg,mult.GI,'-squarek','LineWidth',1.5,'Color','black');
plot(leg,mult.tauc,'-diamondk','LineWidth',1.5,'Color','black');
plot(leg,mult.tauw,'-*','LineWidth',1.5,'Color','black');
%yyaxis right
plot(leg,mult.vg,'-^k','LineWidth',1.5,'Color','red');
plot(leg,mult.MUg,'-ok','LineWidth',1.5,'Color','blue');
xlim([0 40]);
set(gca,'FontName',"Times New Roman",'YColor','black');
set(gca,'FontSize',12);
hold off;

nexttile(tcl);
hold on;
title('Present value unemployment multiplier','Color','black','FontSize',13, 'Interpreter', 'latex');
plot(leg,umult.GC,'-k','LineWidth',1.5,'Color','black');
plot(leg,umult.GI,'-squarek','LineWidth',1.5,'Color','black');
plot(leg,umult.tauc,'-diamondk','LineWidth',1.5,'Color','black');
plot(leg,umult.tauw,'-*','LineWidth',1.5,'Color','black');
%yyaxis right
plot(leg,umult.vg,'-^k','LineWidth',1.5,'Color','red');
plot(leg,umult.MUg,'-ok','LineWidth',1.5,'Color','blue');
xlim([0 40]);
set(gca,'FontName',"Times New Roman",'YColor','black');
set(gca,'FontSize',12);
hold off;
% 将图例置于底部外侧，水平放置
lgd = legend( {'Government consumption',...
             'Public investment', ...
             'Consumption tax cut', ...
             'Labor income tax cut', ...
             'Public vacancy' ...
             'Public wage markup'}, ...
              'Location', 'southoutside', 'Orientation', 'horizontal','NumColumns', 3); 
lgd.FontSize = 12; % 设置图例的字体大小
lgd.FontName = 'Times New Roman'; % 设置图例的字体
fig(2).Position = [0 0 width height]; %[left bottom width height]
hold off;

%% variance decomp
clc,clear; close all;
cd C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_est
load dataset_hp.mat
cd C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_afest

table_list.y_decomp=readtable('est_search_gov_emp_shock_decomposition.xls','Sheet','yobs','VariableNamingRule','preserve');
table_list.n_decomp=readtable('est_search_gov_emp_shock_decomposition.xls','Sheet','nob','VariableNamingRule','preserve');
table_list.pi_decomp=readtable('est_search_gov_emp_shock_decomposition.xls','Sheet','piobs','VariableNamingRule','preserve');
N=length(table_list.y_decomp.Decomposition);
T_list=(datetime(1985,1,1)+calquarters(0:N-1))';
series=[yobs,nob,piobs];

table_list.y_decomp.epstaur=table_list.y_decomp.epstaur/5;
table_list.n_decomp.epstaur=table_list.n_decomp.epstaur/5;
table_list.pi_decomp.epstaur=table_list.pi_decomp.epstaur/5;

table_list.y_decomp.epsgi=table_list.y_decomp.epsgi/1.2;
table_list.n_decomp.epsgi=table_list.n_decomp.epsgi/1.2;
table_list.pi_decomp.epsgi=table_list.pi_decomp.epsgi/1.2;

table_list.y_decomp.epstauc=table_list.y_decomp.epstauc*2;
table_list.n_decomp.epstauc=table_list.n_decomp.epstauc*2;
table_list.pi_decomp.epstauc=table_list.pi_decomp.epstauc*2;

table_list.y_decomp.epset=table_list.y_decomp.epset/1.5;
table_list.n_decomp.epset=table_list.n_decomp.epset/1.5;
table_list.pi_decomp.epset=table_list.pi_decomp.epset/1.5;

table_list.y_decomp.epsnu=table_list.y_decomp.epsnu/2;
table_list.n_decomp.epsnu=table_list.n_decomp.epsnu/2;
table_list.pi_decomp.epsnu=table_list.pi_decomp.epsnu/2;


table_name={'y_decomp','n_decomp','pi_decomp'};
for ii=1:length(table_name)
fig(ii)=figure('Name',append(table_name{ii},'_','figure'));
fig(ii).Position=[100 100 400*1.618 400];
eval(sprintf('table_list.(table_name{%d}).time=T_list;',ii));
eval(sprintf( 'data=table2array(table_list.(table_name{%d})(:,[2:15]));',ii));
yhat=series(:,2);
B=bar(T_list,data(:,[1:14]),'stacked','FaceColor','flat'); hold on;
B(1).FaceColor='r';
B(2).FaceColor='b';
B(3).FaceColor='k';
B(5).FaceColor='m';
B(6).FaceColor='y';
plot(T_list,yhat,'LineWidth',1.5,'Color','red');
xtickformat("yyyy-QQQ")
lg=legend({'Government consumption shock','Public investment shock','Social security shock','Labor income tax shock','Capital income tax shock','Consumption tax shock','Public wage markup shock','Public vacnacy shock','Monetary policy shock','TFP shock','Matiching efficiency shock','Labor seperation shock','Labor supply shock','Price-markup shock'});
lg.Orientation = 'horizontal';
lg.NumColumns=3;
lg.Location='southoutside';
set(gca,'FontSize',11);
set(gca,'FontName',"Times New Roman");
fig(ii).Position = [0 0 450*1.618 450];
hold off;
end


fig(4)=figure('Name',append(table_name{1},'_','fiscal_figure'));
table_list.y_decomp.time=T_list;
yhat=series(:,2);
data_y=table2array(table_list.y_decomp(:,[2:15]));
B=bar(T_list,data_y(:,[1:2 4:8]),'stacked','FaceColor','flat'); hold on;
B(1).FaceColor='r';
B(2).FaceColor='b';
B(3).FaceColor='k';
B(5).FaceColor='m';
B(6).FaceColor='y';
plot(T_list,yhat,'LineWidth',1.5,'Color','red');
xtickformat("yyyy-QQQ")
lg=legend({'Government consumption shock','Public investment shock','Labor income tax shock','Capital income tax shock','Consumption tax shock',"Public wage markup shock","Public vacancy shock"});
lg.Orientation = 'horizontal';
lg.NumColumns=3;
lg.Location='southoutside';
set(gca,'FontSize',11);
set(gca,'FontName',"Times New Roman");
fig(4).Position=[100 100 450*1.618 450];
hold off;


% 提取年份
years = year(table_list.y_decomp.time);
% 初始化年度数据表格
uniqueYears = unique(years);
variable_names=table_list.y_decomp.Properties.VariableNames(1:end);
y_decomp_annual = array2table(zeros(length(uniqueYears), width(table_list.y_decomp)), 'VariableNames',variable_names);
y_decomp_annual.Year = uniqueYears;
%年度加总
variable_names([1 16:end])=[];
for varIdx = 1:length(variable_names)
    varName = variable_names{varIdx};
    for yearIdx = 1:length(uniqueYears)
        y_decomp_annual{yearIdx, varName} = sum(abs(table_list.y_decomp{years==uniqueYears(yearIdx),varName}));
    end
end

for i = 1:height(y_decomp_annual)
    % 计算该年份所有变量的总和
    total = sum(y_decomp_annual{i,variable_names});  % 变量列名列表
    % 计算每个变量的占比
    for varName = variable_names
        y_decomp_annual{i, ['Proportion_' varName{1}]} = y_decomp_annual{i, varName} / total*100;
    end
end

plot_varnames={'Proportion_epsgc','Proportion_epsgi', ...
                'Proportion_epstauw','Proportion_epstaur' ...
                ,'Proportion_epstauc','Proportion_epsvg', ...
                'Proportion_epsMUg'};
proportions = y_decomp_annual{:,plot_varnames};
mean(proportions)
sum(mean(proportions))
fig(5)=figure('Name','relative_importance');
hold on;
bar(y_decomp_annual.Year, proportions, 'stacked');
ylim([0 100])
ytickformat('%.0f%%');
hold off;
lg=legend({'Government consumption (4.34%)','Public investment (3.99%)','Labor income tax (0.54%)','Capital income tax (2.06%)','Consumption tax (4.56%)',"Public vacancy (7.16%)","Public wage markup (1.39%)"});
lg.Orientation = 'horizontal';
lg.NumColumns=2;
lg.Location='southoutside';
set(gca,'FontSize',11);
set(gca,'FontName',"Times New Roman");
fig(5).Position=[100 100 450*1.618 450];
