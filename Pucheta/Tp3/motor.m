% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N� 3 - Caso de estudio 1 - 
%   Inciso 1  
% Implementar un sistema en variables de estado que controle el �ngulo del motor, 
% para consignas de pi/2 y �pi/2 cambiando cada 5 segundos y que el 
% TL de 1,5 aparece s�lo para pi/2, para -pi/2 es nulo. 
% Hallar el valor de integraci�n Euler adecuado. 
% El objetivo es mejorar la din�mica del controlador que muestra la Fig. 1
% verificando el resultado con las curvas del archivo xlsx adjunto
%%
% close all; clear all; clc;
color='r';
% color='m';

%Se importan los datos de la tabla excel
valores=xlsread('Curvas_Medidas_Motor_2023.xls'); 
t_excel=valores(1:end,1);    %tiempo
phi_excel=valores(1:end,2);  %angulo
w_excel=valores(1:end,3);    %velocidad angular
ia_excel=valores(1:end,4);   %corriente
v_excel=valores(1:end,5);    %tension
tl_excel=valores(1:end,6);   %torque

figure(1)
subplot(3,2,1);hold on;
plot(t_excel,phi_excel, 'r' ,'LineWidth',1.5);title('Angulo  \phi [rad]'); grid on;hold on; 
subplot(3,2,2);hold on;
plot(t_excel,w_excel, 'r' ,'LineWidth',1.5);title('Velocidad Angular  \omega [rad/s]');grid on;hold on; 
subplot(3,2,3);hold on;
plot(t_excel,ia_excel, 'r' ,'LineWidth',1.5);title('Corriente  Ia [A]');grid on;hold on; 
subplot(3,2,4);hold on;
plot(t_excel,v_excel, 'r' ,'LineWidth',1.5);title('Tensi�n  V [V]');grid on;hold on;
subplot(3,1,3);hold on;
plot(t_excel,tl_excel, 'r' ,'LineWidth',1.5);title('Torque  TL [N/m]');grid on;hold on;



%Se declaran los par�metros del sistema
LAA = .56; 
J = .0019;
RA = 1.35;
Bm = .000792; % figura 6.9 Kuo (version en ingles) cap 6.2.1
Ki = 1;
Km = 1;

%TIEMPO CONTINUO
%Matrices a lazo abierto
Ac=[-RA/LAA    -Km/LAA    0; %x1 = Corriente
      Ki/J      -Bm/J     0; %x2 = Velocidad angular
       0         1        0] %x3 = Angulo
      
Bc=[1/LAA;
      0;
      0]
      
Cc=[0 0 1;
    1 0 0]                  %Salida theta

Dc=[0]

%val= eig(Ac)

% val=real(1/eig(Ac))

% tR=log(0.95)/val(2); %din�mica mas r�pida
% tR=tR/10;
% tL=log(0.05)/val(3); %din�mica mas lenta
% tL=tL*10;


%Tiempos
Ts = 3e-3;         %Tiempo de muestreo
T = 20;              %Tiempo de simulaci�n � M�s lento
At = Ts;            %Tiempo de integraci�n � M�s r�pido
Kmax = (T/At);

%Definici�n de variables
titaRef = (-pi/2);   %Valor de referencia inicial
Tl =1.5; %0.1035 %este valor es el promedio del excel        %Torque (solo para pi/2)
tc = 5;        %Tiempo en el que se pasa la referencia de pi/2 a -pi/2
est = 0;

%Asignaciones
t = 0:At:T;
ia = 0:At:T;
w = 0:At:T;
w_p = 0:At:T;
theta = 0:At:T;

u = linspace(0,0,Kmax+1);
ref = linspace(0,0,Kmax+1);
torque = linspace(0,0,Kmax+1);
%Condiciones Iniciales
ia(1) = 0;
w(1) = 0;
w_p(1) = 0;
theta(1) = 0;
TL = 0;        %El torque comienza en Tl, ya que la referencia inicial es -pi/2
u(1) = 0;

x = [0,0,0];      %[ia, tita, omega] condiciones inicales
x_p = [0,0,0];    %[ia_p, tita_p, omega_p]
x_hat = [0;0;0];


%Pasaje de tiempo continuo a tiempo discreto
sys_c = ss(Ac,Bc,Cc,Dc);
sys_d = c2d(sys_c,Ts,'zoh');

%TIEMPO DISCRETO
%Matrices de lazo abierto 
A = sys_d.a;
B = sys_d.b;
C = sys_d.c;
 
%Matriz de Ma y Mc
% Ma = [B A*B A^2*B A^3*B]; Alcanzabilidad = rank(Ma)    %Alcanzabilidad
% Mc = [B A*B A^2*B A^3*B A^4]; Controlabilidad = rank(Mc) %Controlabilidad

% Para el integrador
Cref = Cc(1,:);
AA=[A,zeros(3,1);-Cref*A,eye(1)]
BB=[B;-Cref*B]


%Par�metros del controlador LQR
% d = [1e-5 1e-6 1.9e2 0.001];%Corriente, �ngulo, Velocidad angular
d = [1 1 1 0.1];%Corriente, velocidad , �ngulo
Q = diag(d);
R = 1e-2;            %Dimensionar para que no pase de 24 v
Kk = dlqr(AA, BB, Q, R);
K= Kk(1:3);
Ki = -Kk(4);


%Punto N�2: se implementa un controlador con observador
%Observador
Ao = A';
Bo = Cc';
Co = B';
%Par�metros del controlador LQR observador
do = [1e1 10 .09e0];%Corriente, �ngulo, Velocidad angular
Qo = diag(do); 
Ro=[1e1    0  ;
      0    1e1];
Ko = (dlqr(Ao,Bo,Qo,Ro))'

%Ganancia de prealimentaci�n
Gj = inv(Cc(1,:)*inv(eye(3)-A+B*K)*B);

%Simulaci�n
ii = 1;
kk = 0;
% psi(1)=0;
v(1) = 0;
while(ii<(Kmax))
    kk=kk+At;
    if(kk>tc)
        titaRef=titaRef*(-1);
        if(est==0)
            TL=Tl;
            est=1;
        else
            TL=0;
            est=0;
        end
        kk=0;
    end
    ref(ii)=titaRef;
    torque(ii)=TL ;
    estado=[ia(ii); w(ii); theta(ii)];
    
    y_sal=Cc*estado;
    v(ii+1)=v(ii)+ref(ii)-y_sal(1);
    
    %Ley de control
%     u(ii)=-K*estado+Ki*v(ii+1); %Sin observador  %+Gj*ref(ii)
    u(ii) = -K*x_hat+Ki*v(ii+1); %Con observador
    
    %Punto N�4: se aumenta la zona muerta hasta alcanzar valores inaceptables
    zona_muerta=1;
    
    %Zona Muerta
    if(abs(u(ii))<zona_muerta)
        u(ii)=0;
    else
        u(ii)=sign(u(ii))*(abs(u(ii))-zona_muerta);
    end
    
    %----------------------------------------------
    
%     y_sal = Cc*estado;
    y_sal_o=Cc*x_hat;
    
    %Integracion de Euler
    ia_p = -(RA/LAA)*ia(ii) - (Km/LAA)*w(ii) + (1/LAA)*u(ii);
    w_p = (Ki/J)*ia(ii) -(Bm/J)*w(ii) + (1/J)*TL;
    tita_p = w(ii);
    x_hat = A*x_hat+B*u(ii)+Ko*(y_sal - y_sal_o);
    ia(ii+1) = ia(ii)+At*ia_p;
    w(ii+1) = w(ii)+At*w_p;
    theta(ii+1) = theta(ii)+At*w(ii);
    ii=ii+1;
  
end


%Gr�ficos

figure(2)
subplot(3,2,1);hold on;
plot(t,theta, color ,'LineWidth',1.5);title('Angulo  \phi [rad]'); grid on;hold on;
plot(t,ref,'g','LineWidth',1.5)
subplot(3,2,2);hold on;
plot(t,w, color ,'LineWidth',1.5);title('Velocidad Angular  \omega [rad/s]');grid on;hold on; 
subplot(3,2,3);hold on;
plot(t,ia, color ,'LineWidth',1.5);title('Corriente  Ia [A]');grid on;hold on; 
subplot(3,2,4);hold on;
plot(t,u, color ,'LineWidth',1.5);title('Acci�n de control  V [V]');grid on;hold on;
subplot(3,1,3);hold on;
plot(t,torque, color ,'LineWidth',1.5);title('Torque  TL [N/m]');grid on;hold on;

figure(3)
plot(theta,w, color ,'LineWidth',1.5);title('Plano de Fases'); grid on;hold on;