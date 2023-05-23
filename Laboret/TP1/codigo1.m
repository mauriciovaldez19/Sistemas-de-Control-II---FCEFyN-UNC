clear all
close all
clc

% Datos: Polo1=-2       Polo2=0     Cero=-10
%Ganancia=10    Sobrepaso=5    tiempo2%=3   error=0 periodo=0.15
%% Definicion FdT continua %%
G=zpk([-10],[-2 0],[10]);

%%  Definicion FdT discreta
Tm=0.15;
Gd=c2d(G,Tm,'zoh');
%% Mapa de Polos y Ceros
pzmap(G) 
pzmap(Gd)

%% Aumento x 10 el periodo de muestreo
Gd1=c2d(G,10*Tm,'zoh');
pzmap(Gd1)

%% Entrada escalon a G y Gd
step(G)
step(Gd)

%Analisis sistema discreto realimentado
Kp=dcgain(Gd)
error=1/Kp
F=feedback(Gd,1)
step(F)

t=0:Tm:100*Tm %genera rampa
lsim(F,t,t)


%% Analisis a lazo cerrado con realimentacion unitaria
G=zpk([-10],[-2 0],[10]);
Tm=0.15
Gd=c2d(G,Tm,'zoh')
Gd1=c2d(G,10*Tm,'zoh')
rlocus(G)
rlocus(Gd)

%Aumentando 10 veces el tiempo de muestreo original
rlocus(Gd1)
