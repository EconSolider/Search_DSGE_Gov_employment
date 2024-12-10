%% price stick
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


omegaw_list=linspace(0.1,0.9,30);
omegaw=0.5956;

for ii=1:length(omegaw_list)
omegap=omegaw_list(ii);
dynare afest_search_gov_emp nolog
wreve=oo_.irfs.tauw_epstauw* ...
      (steady_state.wp*steady_state.np*steady_state.hp ...
        +steady_state.wg*steady_state.ng*steady_state.hg);
creve=oo_.irfs.tauc_epstauc*steady_state.C;
    mult.vg_1(ii)=abs(multi(oo_.irfs.Y_epsvg(1:end),oo_.irfs.deltb_epsvg(3:end),1,0.97));
    mult.MUg_1(ii)=abs(multi(oo_.irfs.Y_epsMUg(1:end),oo_.irfs.deltb_epsMUg(3:end),1,0.97));
    mult.GC_1(ii)=abs(multi(oo_.irfs.Y_epsgc(1:end),oo_.irfs.deltb_epsgc(2:end),1,0.97));
    mult.GI_1(ii)=abs(multi(oo_.irfs.Y_epsgi(1:end),oo_.irfs.gi_epsgi(2:end),1,0.97));
    mult.tauc_1(ii)=abs(multi(oo_.irfs.Y_epstauc(1:end),creve(2:end),1,0.97));
    mult.tauw_1(ii)=abs(multi(oo_.irfs.Y_epstauw(1:end),wreve(2:end),1,0.97))/5;
    
    mult.vg_40(ii)=abs(multi(oo_.irfs.Y_epsvg(1:end),oo_.irfs.deltb_epsvg(3:end),40,0.97));
    mult.MUg_40(ii)=abs(multi(oo_.irfs.Y_epsMUg(1:end),oo_.irfs.deltb_epsMUg(3:end),40,0.97));
    mult.GC_40(ii)=abs(multi(oo_.irfs.Y_epsgc(1:end),oo_.irfs.deltb_epsgc(2:end),40,0.97));
    mult.GI_40(ii)=abs(multi(oo_.irfs.Y_epsgi(1:end),oo_.irfs.gi_epsgi(2:end),40,0.97));
    mult.tauc_40(ii)=abs(multi(oo_.irfs.Y_epstauc(1:end),creve(2:end),40,0.97));
    mult.tauw_40(ii)=abs(multi(oo_.irfs.Y_epstauw(1:end),wreve(2:end),40,0.97))/5;

    umult.vg_1(ii)=-abs(multi(oo_.irfs.n_epsvg(2:end),oo_.irfs.deltb_epsvg(3:end),1,0.97))/2.5;
    umult.MUg_1(ii)=-multi(oo_.irfs.n_epsMUg(2:end),oo_.irfs.deltb_epsMUg(3:end),1,0.97)*30;
    umult.GC_1(ii)=-abs(multi(oo_.irfs.n_epsgc(2:end),oo_.irfs.deltb_epsgc(2:end),1,0.97))*30;
    umult.GI_1(ii)=-abs(multi(oo_.irfs.n_epsgi(2:end),oo_.irfs.gi_epsgi(2:end),1,0.97))*30;
    umult.tauc_1(ii)=-multi(oo_.irfs.n_epstauc(2:end),creve(2:end),1,0.97)*30;
    umult.tauw_1(ii)=-multi(oo_.irfs.n_epstauw(2:end),wreve(2:end),1,0.97)/3;
    
    umult.vg_40(ii)=-abs(multi(oo_.irfs.n_epsvg(2:end),oo_.irfs.deltb_epsvg(3:end),40,0.97))/2.5;
    umult.MUg_40(ii)=-multi(oo_.irfs.n_epsMUg(2:end),oo_.irfs.deltb_epsMUg(3:end),40,0.97)*30;
    umult.GC_40(ii)=-abs(multi(oo_.irfs.n_epsgc(2:end),oo_.irfs.deltb_epsgc(2:end),40,0.97))*30;
    umult.GI_40(ii)=-abs(multi(oo_.irfs.n_epsgi(2:end),oo_.irfs.gi_epsgi(2:end),40,0.97))*30;
    umult.tauc_40(ii)=-multi(oo_.irfs.n_epstauc(2:end),creve(2:end),40,0.97)*30;
    umult.tauw_40(ii)=-multi(oo_.irfs.n_epstauw(2:end),wreve(2:end),40,0.97)/3;
end

%% %% relation between multiplier and omegaw
width=600*1.618;
height=500;
fig(1)=figure('Name','relation between multiplier and $$\omega_w$$');
tcl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile(tcl);
hold on;
plot(omegaw_list,mult.vg_1,'-^k','LineWidth',1.5,'Color','red');
plot(omegaw_list,mult.MUg_1,'-ok','LineWidth',1.5,'Color','blue');
plot(omegaw_list,mult.GC_1,'-k','LineWidth',1.5,'Color','black');
plot(omegaw_list,mult.GI_1,'-squarek','LineWidth',1.5,'Color','black');
plot(omegaw_list,mult.tauc_1,'-diamondk','LineWidth',1.5,'Color','black');
plot(omegaw_list,mult.tauw_1,'-*','LineWidth',1.5,'Color','black');
hold off;
ylim([0 1.5]);
xline(0.6986,'--','color','red');
title('PV(Y_{0})','Color','black','FontSize',12);
xlabel('Price rigidity parameter \omega^p');
set(gca,'FontSize',12);
set(gca,'FontName',"Times New Roman");

nexttile(tcl);
hold on;
plot(omegaw_list,mult.vg_40,'-^k','LineWidth',1.5,'Color','red');
plot(omegaw_list,mult.MUg_40,'-ok','LineWidth',1.5,'Color','blue');
plot(omegaw_list,mult.GC_40,'-k','LineWidth',1.5,'Color','black');
plot(omegaw_list,mult.GI_40,'-squarek','LineWidth',1.5,'Color','black');
plot(omegaw_list,mult.tauc_40,'-diamondk','LineWidth',1.5,'Color','black');
plot(omegaw_list,mult.tauw_40,'-*','LineWidth',1.5,'Color','black');
hold off;
ylim([0 1.5]);
xline(0.6986,'--','color','red');
title('PV(Y_{40})','Color','black','FontSize',12);
xlabel('Price rigidity parameter \omega^p');
set(gca,'FontSize',12);
set(gca,'FontName',"Times New Roman");

nexttile(tcl);
hold on;
plot(omegaw_list,umult.vg_1,'-^k','LineWidth',1.5,'Color','red');
plot(omegaw_list,umult.MUg_1,'-ok','LineWidth',1.5,'Color','blue');
plot(omegaw_list,umult.GC_1,'-k','LineWidth',1.5,'Color','black');
plot(omegaw_list,umult.GI_1,'-squarek','LineWidth',1.5,'Color','black');
plot(omegaw_list,umult.tauc_1,'-diamondk','LineWidth',1.5,'Color','black');
plot(omegaw_list,umult.tauw_1,'-*','LineWidth',1.5,'Color','black');
hold off;
ylim([-0.5 0.2]);
xline(0.6986,'--','color','red');
title('PV(u_{1})','Color','black','FontSize',12);
xlabel('Price rigidity parameter \omega^p');
set(gca,'FontSize',12);
set(gca,'FontName',"Times New Roman");

nexttile(tcl);
hold on;
plot(omegaw_list,umult.vg_40,'-^k','LineWidth',1.5,'Color','red');
plot(omegaw_list,umult.MUg_40,'-ok','LineWidth',1.5,'Color','blue');
plot(omegaw_list,umult.GC_40,'-k','LineWidth',1.5,'Color','black');
plot(omegaw_list,umult.GI_40,'-squarek','LineWidth',1.5,'Color','black');
plot(omegaw_list,umult.tauc_40,'-diamondk','LineWidth',1.5,'Color','black');
plot(omegaw_list,umult.tauw_40,'-*','LineWidth',1.5,'Color','black');
ylim([-0.5 0.2]);
xline(0.6986,'--','color','red');
title('PV(u_{40})','Color','black','FontSize',12);
xlabel('Price rigidity parameter \omega^p');
set(gca,'FontSize',12);
set(gca,'FontName',"Times New Roman");
hold off;
fig(1).Position = [0 0 width height];
lgd = legend( {'Public vacancy', ...
             'Public wage markup', ...
             'Government consumption',...
            'Public investment', ...
             'Consumption tax cut', ...
             'Labor income tax cut'});
lgd.Orientation='horizontal';
lgd.NumColumns=3;
lgd.Layout.Tile = 'south';
lgd.FontSize = 12; % 设置图例的字体大小
lgd.FontName = 'Times New Roman'; % 设置图例的字体
