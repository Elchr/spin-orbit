clc; clear;

syms eccm e m a q w
eccm=m+e*sin(m+e*sin(m+e*sin(eccm)));
ecaser=series(eccm,e,'Order',3); 
rser=combine(series(a*(1-e*cos(ecaser)),e,'Order',3),'sincos');
cosfser=simplify(combine(series((cos(eccm)-e)/(1-e*cos(eccm)),e,'Order',3),'sincos'));
sinfser=(simplify(rewrite(simplify(series(sqrt(1-e^2)*sin(eccm)/(1-e*cos(eccm)),e,'Order',3)),'sin')));

total=simplify(1/rser^3*(sin(2*(w/2+m))*(cosfser^2-sinfser^2)-cos(2*(w/2+m))*2*sinfser*cosfser));
totalser=simplify(collect(simplify(series(total,e,'Order',3)),e));

