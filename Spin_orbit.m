clc; clear;
%syms m w e a 
%totalser=(-7*e*sin(m-w)-17*e^2*sin(2*m-w)+2*sin(w)-5*e^2*sin(w)-e*sin(m+w))/(2*a^3);

M=0.0002858; al=0.00334; G=4*(pi)^2;
n=1000;
heta=sqrt(G*M/al^3); T=(2*pi)/sqrt(G*M/al^3);
dt=T/1000;
total=n*T;
nstep=fix(total/dt);

m2=zeros(nstep,21);
totalserval=zeros(nstep,21);
wold=zeros(nstep,21);
pwold=zeros(nstep,21);
wnew=zeros(nstep,21);
pwnew=zeros(nstep,21);
t=zeros(nstep,21);
k=1;
eps=0.7;
el=0.8;

for h=-10:1:10
    wold(1,k)=0;
    pwold(1,k)=h*3*heta*eps/10;
for i=1:1:nstep;
        
    t(i,k)=i*dt;
    m2(i,k)=2*pi*t(i,k)/T;

  totalserval(i,k)=-7*el*sin(m2(i,k)-wold(i,k))/(2*al^3)-17*el^2*sin(2*m2(i,k)-wold(i,k))/(2*al^3)+sin(wold(i,k))/al^3-5*el^2*sin(wold(i,k))/(2*al^3)-el*sin(m2(i,k)+wold(i,k))/(2*al^3) ; 
  pwnew(i,k)=+pwold(i,k)-G*M*eps^2*totalserval(i,k)*dt;%*wold(i,k);
  wnew(i,k)=mod((wold(i,k)+pwnew(i,k)*dt)+pi,2*pi)-pi;
       
  wold(i+1,k)=wnew(i,k);
  pwold(i+1,k)=pwnew(i,k);
end
k=k+1;
end 

b1x=zeros(1000,21);
b1y=zeros(1000,21);

for jj=1:1:21
for j=1:1:n
b1x(j,jj)=wold(j*(T/dt),jj);
b1y(j,jj)=pwold(j*(T/dt),jj); 
end
end

%figure(2)
%for i=1:1:21
%plot (b1x(:,i),b1y(:,i),'.k')
%hold on

figure(1)
plot (b1x,b1y,'.k')
xlabel('\Psi')
ylabel('$\dot{\Psi}$','interpreter','latex')
grid on
title('Διάγραμμα Φάσης')
legend('e=0.5','ε=0.7')



