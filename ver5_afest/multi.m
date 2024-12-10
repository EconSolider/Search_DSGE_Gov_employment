function PV_multi=multi(Y,G,x,beta)
% input 
% Y:IRF of GDP
% G:IRF of government expenditure
% beta: discount rate
% x: periods after policy shock
discount_Y=zeros(x,1);
discount_G=zeros(x,1);
for ii=1:x
discount_Y(ii)=beta^(ii-1) * Y(ii);
discount_G(ii)=beta^(ii-1) * G(ii);
end
PV_Y=sum(discount_Y);
PV_G=sum(discount_G);
PV_multi=PV_Y/PV_G;
end