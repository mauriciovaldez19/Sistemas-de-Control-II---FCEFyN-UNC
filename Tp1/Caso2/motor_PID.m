% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 1 - Caso de estudio 2 - Inciso 4  
%   Implementar un PID en tiempo discreto para que el ángulo del motor permanezca en una 
%   referencia de 1radian. (Tip: partir de KP=0,1; Ki=0,01; KD=5).

%%
clc;clear;close all;
X=-[0; 0; 0; 0];ii=0;t_etapa=1e-7;tF=.6;wRef=1;
color_='m';

%Constantes del PID
Kp=0.5;Ki=8;Kd=0; %definitivas
%Kp=1;Ki=8;Kd=1; %de puebas
color_='r';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(tF/t_etapa,1);u=0;

for t=0:t_etapa:tF
ii=ii+1;k=ii+2;
Tl=2.13e-5/100; % Tl se busca un valor mucho mas chico 
                % que el limite encontrado en inciso 1
X=modmotorpunto4(t_etapa, X, u,Tl);

e(k)=wRef-X(4); %ERROR
u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
    if u>12         %limito accion de control a +-12
        u=12;
    end
    if u<-12
        u=-12;
    end

x1(ii)=X(1);%Omega
x2(ii)=X(2);%wp
x3(ii)=X(3);%ia
x4(ii)=X(4);%tita
acc(ii)=u; %accion de control
end
t=0:t_etapa:tF;


figure(1)
subplot(4,1,1);hold on;
plot(t,x4,color_);title('Angulo tita rad');hold on;

subplot(4,1,2);hold on;
plot(t,acc,color_);title('Accion de control');hold on;

subplot(4,1,3);hold on;
plot(t,x3,color_);title('Ia');hold on;

subplot(4,1,4);hold on;
plot(t,x1,color_);title('\omega_t');hold on;


%motor
function [X]=modmotorpunto4(t_etapa, xant, accion,Tl)
Laa=366e-6; J=5e-9;Ra=55.6;B=0;Ki=6.49e-3;Km=6.53e-3;Va=accion;
h=1e-7;
TL=Tl;   
omega= xant(1);     
wp= xant(2);        
ia=xant(3);
tita = xant(4);
for ii=1:t_etapa/h
wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
iap=(-Ra*ia-Km*omega+Va)/Laa;
wp=wp+h*wpp;
wp=wp-(TL/J);
ia=ia+iap*h;
omega = omega + h*wp;
tita = tita + h*omega;
end
X=[omega,wp,ia,tita];
end