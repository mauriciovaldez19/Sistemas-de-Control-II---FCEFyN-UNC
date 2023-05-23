% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Ing. Pucheta
% Alumno: Valdez Benavidez, Mauricio Luciano
% Consigna 2 - Tp N° 1
%%
clc;clear;close all;
X=-[0; 0; 0; 0];ii=0;t_etapa=1e-7;tF=.6;wRef=1;
color_='m';

%Constantes del PID
Kp=1;Ki=0;Kd=0;color_='r';

Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(tF/t_etapa,1);u=0;

for t=0:t_etapa:tF
ii=ii+1;k=ii+2;
Tl=0;
     if(t>=0.1501)  % tiempo tabla excel
        Tl=4.5e-2; % Tl que se ajusta a la rta de los datos del excel
     end
     
X=modmotorpunto2(t_etapa, X, u,Tl);
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

% %Cargo valores para la Comparacion
% valores=xlsread('Curvas_Medidas_Motor_2023.xlsx'); 
% tt=valores(1:end,1);
% W=valores(1:end,2);


figure(1)
subplot(4,1,1);hold on;
plot(t,x4,color_);title('Angulo tita rad');hold on;

subplot(4,1,2);hold on;
plot(t,acc,color_);title('Accion de control');hold on;

subplot(4,1,3);hold on;
plot(t,x3,color_);title('Ia');hold on;

subplot(4,1,4);hold on;
plot(t,x1,color_);title('Salida y, \omega_t');hold on;



%motor
function [X]=modmotorpunto2(t_etapa, xant, accion,Tl)
Laa=506.94e-6; J=3.963e-6;Ra=99.7224;B=0;Ki=16.52;Km=(1/16.52);
Va=accion;
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