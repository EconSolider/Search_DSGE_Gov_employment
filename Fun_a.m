function F=Fun_a(x,alph,pm,y,hp,mu,bet,I,gi,gc,u,z_wp,vp,vg,Y,b,i,T,ng,tauw,MUg,np,hg,taur,r,K,tauc)
C=x(1);
lp=x(2);
wp=x(3);
lg=0.15*lp;

% F(1)=(1-alph)*pm*y-wp*hp+(1-chi)*lp/mu ...
%     -lp/bet/mu;
F(1)=(1-alph)*pm*y-wp*hp ...
     -lp/bet/mu;

F(2)=C+I+gi+gc+u*z_wp*wp+vp*lp+vg*lg ...
    -Y;
F(3)=i*b+u*z_wp*wp+gi+gc+(1+MUg)*wp*ng*hg+vg*lg-T ...
    -tauw*wp*np*hp-tauw*(1+MUg)*wp*ng*hg-taur*r*K-tauc*C ...
    -b;
end