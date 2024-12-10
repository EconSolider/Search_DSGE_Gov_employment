%% price stickness and policy efficiency check
clc,clear;
close all;
cd C:\Users\lzkzz\Desktop\matlearn\Dynare_train\search_gov_emp\ver5_afest

omegap=0.6986;
omegaw=0.5956;
dynare afest_search_gov_emp nolog
clearvars -except oo_ M_
for j=1:length(oo_.var_list)
   eval(sprintf('steady_state.(oo_.var_list{%d})=oo_.steady_state(%d);',j,j));
end

omegap_list=[0.7 0.9];
omegaw=0.5956;
for ii=1:length(omegap_list)
omegap=omegap_list(ii);
dynare afest_search_gov_emp nolog
oo_list(ii)=oo_;
end

oo_a=oo_list(1);
oo_b=oo_list(2);

for j=1:length(oo_.var_list)
   eval(sprintf('steady_state.(oo_.var_list{%d})=oo_.steady_state(%d);',j,j));
end
irfnames=fieldnames(oo_.irfs);
split_irfnames=split(irfnames,'_');
varlist=split_irfnames(:,1);
pctvar=["pi" "pistar" "r" "i" "np" "ng" "n" "u" "fp" "fg" ... 
    "thet" "taur" "tauc" "tauw" "chi" "deltb" "MUg" "vg" "vp"];
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
     eval(sprintf('irfs_a.(irfnames{%d})=oo_a.irfs.(irfnames{%d})(2:end)/abs(steady_state.(varlist{%d}))*100;',j,j,j));
     eval(sprintf('irfs_b.(irfnames{%d})=oo_b.irfs.(irfnames{%d})(2:end)/abs(steady_state.(varlist{%d}))*100;',j,j,j));
 else
     eval(sprintf('irfs_a.(irfnames{%d})=oo_a.irfs.(irfnames{%d})(2:end)*100;',j,j));
     eval(sprintf('irfs_b.(irfnames{%d})=oo_b.irfs.(irfnames{%d})(2:end)*100;',j,j));
 end
end

irfs_a.C_epsgc=-oo_a.irfs.C_epsgc(3:end);
irfs_b.C_epsgc=-oo_b.irfs.C_epsgc(3:end);

irfs_a.n_epstauc=-oo_a.irfs.n_epstauc(3:end);
irfs_a.np_epstauc=-oo_a.irfs.np_epstauc(3:end);
irfs_b.n_epstauc=-oo_b.irfs.n_epstauc(3:end);
irfs_b.np_epstauc=-oo_b.irfs.np_epstauc(3:end);

irfs_a.n_epsgi=-oo_a.irfs.n_epsgi(3:end);
irfs_a.np_epsgi=-oo_a.irfs.np_epsgi(3:end);
irfs_b.n_epsgi=-oo_b.irfs.n_epsgi(3:end);
irfs_b.np_epsgi=-oo_b.irfs.np_epsgi(3:end);

irfs_a.n_epstauw=-oo_a.irfs.n_epstauw(3:end);
irfs_a.np_epstauw=-oo_a.irfs.np_epstauw(3:end);
irfs_a.wlower_epstauw=-oo_a.irfs.wlower_epstauw(3:end);
irfs_b.n_epstauw=-oo_b.irfs.n_epstauw(3:end);
irfs_b.np_epstauw=-oo_b.irfs.np_epstauw(3:end);
irfs_b.wlower_epstauw=-oo_b.irfs.wlower_epstauw(3:end);

irfs_a.wp_epsvg=-oo_a.irfs.wp_epsvg;
irfs_b.wp_epsvg=-oo_b.irfs.wp_epsvg;
clearvars -except irfs_a irfs_b
%画像設定
width=600*1.618;
height=600;
%% irfs
endo_var_list=["vg" "Y" "C" "I" "np" "hp" "ng" "hg" "n"  ...
               "thet" "fp" "fg" "wp" ...
               "wupper" "wlower" "r" "i" "pi" "pm" ...
               "benefit" "gs" "b" "wreve"];
endo_var_name=["$$v^g$$" "$$Y$$" "$$C$$" "$$I$$" "$$n^p$$" "$$h^p$$" "$$n^g$$" "$$h^g$$" "$$n$$" ...                
               "$$\theta$$" "$$f^p$$" "$$f^g$$" "$$w^p$$" ...
               "$$\overline{w}$$" "$$\underline{w}$$" "$$r$$" "$$i$$" "$$\pi$$" "$$p_m$$"...
               "$$z \times u$$" "$$g^s$$" "$$b$$" "wreve"];
shock_var_list=["epstauc"];


irf_list=[];
for i=1:length(shock_var_list)
    for j=1:length(endo_var_list)
        irf_name=append(endo_var_list(j),'_',shock_var_list(i));
        irf_list=[irf_list; irf_name];
    end
end

fig(1)=figure('Name','government vacancy shock');
tcl = tiledlayout(6,4,'TileSpacing','compact','Padding','compact');
for i=1:length(endo_var_list)
    irf_a=append('irfs_a.',irf_list(i));
    irf_b=append('irfs_b.',irf_list(i));
    nexttile(tcl);
   % x=0:length(eval(irf))-1;
    plot(eval(irf_a),'Color','black','LineWidth',1.5); hold on;
    plot(eval(irf_b),'Color','blue','LineWidth',1.5); hold on;
    yline(0,'--','Color','red','LineWidth',1);
    title(endo_var_name(i),'Color','blue','FontSize',13, 'Interpreter', 'latex');
    set(gca,'FontName',"Times New Roman");
    set(gca,'FontSize',11);
end
fig(1).Position = [0 0 width height]; %[left bottom width height]
hold off;