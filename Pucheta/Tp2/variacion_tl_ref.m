clear all, close all, clc
Tf=12;h=1e-4; pasos=Tf/h; t=0:h:Tf;
t=linspace(0,Tf,pasos);
ref=linspace(0,0,pasos);
tl=linspace(0,0,pasos);
Tref=pi/2;
tc=2;
ii=0;
Tl=1.15e-3;

for i=1:pasos-1
 ii=ii+h;
 if (ii>=tc)
 ii=0;
 Tref=Tref*-1;
 end
 ref(i)=Tref;
 if Tref>=0
 tl(i)=Tl;
 else
     tl(i)=0;
 end
end
figure(1)
plot(t,ref)
figure(2)
plot(t,tl)
