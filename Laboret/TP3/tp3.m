%%Tp3 - Ing Laboret
%Alumno: Valdez Benavidez, Mauricio Luciano

% Datos: Polo1=2   Polo2=-2  Ganancia=10    tiempo2%=4
%%
p1=2;
p2=-2;
K=10;
G=zpk([],[p1 p2],K)
sisotool(G)
% Exportar el controlador (C)

%Kc=0.2978 % ganancia
Kc=0.5
a=2 % cero del controlador con signo invertido
M=1 % ganancia rele 
T=0.1 % histeresis
lineal=1 % simula control lineal
sim('bang_bang_hist_DI_PD')
figure(1)
subplot(2,2,1);plot(tout,yout(:,1),'r','LineWidth',2);grid on;hold on; title('Error'); % error
subplot(2,2,2);plot(yout(:,1),yout(:,3),'b','LineWidth',2);grid on; hold on; title('Plano de fases');
xlabel('Error');ylabel('Derivada del error');% plano de fases: eje x error, eje y derivada del error
subplot(2,1,2);plot(tout,yout(:,2),'m','LineWidth',2);grid on; hold on; title('Acción de control'); % señal de control

M=Kc % ganancia rele = +-ganancia Kc
T=K*Kc/100 % 100 veces menos que la ganancia total
lineal=0 % simula no lineal


T=K*Kc/25 % 25 veces menos que la ganancia total
T=K*Kc/10 % 10 veces menos que la ganancia total

T=K*Kc 