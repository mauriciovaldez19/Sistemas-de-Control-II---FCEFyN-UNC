% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N� 2 - Caso de estudio 1 - 
%   Inciso 1  
% Implementar un sistema en variables de estado que controle el �ngulo del motor, 
% para consignas de pi/2 y �pi/2 cambiando cada 2 segundos y que el 
% TL de 1,15 10-3 aparece s�lo para pi/2, para -pi/2 es nulo. 
% Hallar el valor de integraci�n Euler adecuado. 
% El objetivo es mejorar la din�mica del controlador que muestra la Fig. 1 
%%
clc;clear all;close all;
%Parametros de motor
Laa=5e-3;J=0.004;
Ra=0.2;Bm=0.005;
Ki=6.5e-5;
Km=0.055;

%Matrices de estado
A = [-Ra/Laa -Km/Laa 0;Ki/J -Bm/J 0;0 1 0];
B = [1/Laa 0;0 -1/J;0 0]; %considerando el torque
C = [0 0 1];             %salida posicion

%Matrices ampliadas
Aa=[A zeros(3,1); -C 0];
Ba=[B(:,1); 0];
Ca=[C 0];

%LQR
Q=diag([10 1/100 1/100 10000/2.5]);    R=10;
%Q=diag([60 .001 20 5000]);    R=10;

Kamp=lqr(Aa,Ba,Q,R);

Tf=100;h=1e-4;pasos=Tf/h; t=linspace(0,Tf,pasos);
Ref=linspace(0,0,pasos);
TL=linspace(0,0,pasos);
Tref=pi/2;
Tl=1.15e-3;
tc=20;
ii=0;

for i=1:pasos-1
 ii=ii+h;
 if (ii>=tc)
 ii=0;
 Tref=Tref*-1;
 end
 Ref(i)=Tref;
 if Tref>=0
 TL(i)=Tl;
 else
     TL(i)=0;
 end
end


%condiciones iniciales
ia(1)=0;                %Corriente de armadura x1
theta(1)=0;             %Valocidad angular x2
omega(1)=0;             %Posicion angular x3

estados=[ia(1);
        omega(1);
        theta(1)];

    
Xop=[0 0 0]';
x=[ia(1) omega(1) theta(1)]';

psi(1)=0;
integracion(1)=psi(1);

for i=1:round(Tf/h)
    
    psi_p=Ref(i)-C*estados;
    psi(i)=integracion+psi_p*h;
    
    u(i) = -Kamp(1:3)*estados-Kamp(4)*psi(i);

     %saturador en +-12V
 if u(i)>12
 u(i)=12;
 end
  if u(i)<-12
 u(i)=-12;
 end 
    
    %Variables del sistema lineal
    ia(i)= x(1);
    omega(i)= x(2);
    theta(i)= x(3);
    
    x1_p=-Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
    x2_p=Ki*x(1)/J-Bm*x(2)/J-TL(i)/J;
    x3_p=x(2);
    
    xp=[x1_p; x2_p; x3_p];
    x=x+h*xp;
    
    estados=[ia(i);omega(i);theta(i)];
    integracion=psi(i);
end

figure(2)
subplot(2, 2, 1);
plot(t,ia);
title('Corriente de armadura i_a');
xlabel('Tiempo (seg.)');
ylabel('Corriente (A)');
grid on;

subplot(2, 2, 2);
plot(t,omega);
title('Velocidad angular \omega_r');
xlabel('Tiempo (seg.)');
ylabel('Velocidad angular (rad/s)');
grid on;

subplot(2, 2, 3);
hold on
plot(t,theta);
plot(t,Ref);
hold off
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(2, 2, 4);
hold on
plot(t,u);
hold off
title('Accion de control u_t');
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;