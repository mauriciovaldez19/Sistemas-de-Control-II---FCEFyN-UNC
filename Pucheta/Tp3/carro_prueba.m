%Caso de estudio 3. Sistema no lineal de cuatro variables de estado
close all;clear all; clc;
color = 'r';

%Se declaran los parámetros del sistema
m = 0.1;
Fricc = 0.1; 
l = 1.6; 
g = 9.8;
M = 1.5;

%TIEMPO CONTINUO
%Matrices de lazo abierto en el equilibrio estable
Ac=[0       1               0       0;  %x1=delta - desplazamiento
    0    -Fricc/M         -m*g/M    0;  %x2=delta_p
    0       0               0       1;  %x3=phi - angulo
    0  -Fricc/(l*M)  -g*(m+M)/(l*M) 0];  %x4=phi_p
    
Bc=[  0; 
     1/M;
      0;
    1/(l*M)];
    
Cc=[1 0 0 0;
    0 0 1 0];
   
Dc=[0];


%Condiciones iniciales
alfa(1) = pi;               %Ángulo inicial
ref = 10;               %Posición de referencia           
flag = 0;


%Tiempos
Ts = 1e-2;          
T = 15;            
At = 1e-4;          
Kmax = T/Ts;

%Pasaje de tiempo continuo a tiempo discreto
sys_c = ss(Ac, Bc, Cc, Dc);
sys_d=c2d(sys_c,Ts,'zoh');

%TIEMPO DISCRETO
%Matrices de lazo abierto en el equilibrio estable
A = sys_d.a;
B = sys_d.b;
C = sys_d.c;

%Integrador se mide desplazamiento
Cref = Cc(1,:);

AA=[A,zeros(4,1);-Cref*A,eye(1)];
BB=[B;-Cref*B];

% %Matriz de Ma y Mc
% Ma = [BB AA*BB AA^2*BB AA^3*BB AA^4*BB AA^5*BB]; Alcanzabilidad = rank(Ma);
% Mc = [BB AA*BB AA^2*BB AA^3*BB AA^4*BB AA^5*BB AA^6]; Controlabilidad = rank(Mc); 


%parametros DLQR   de 0 a 10 m
dd = [.1 .1 .1 .1 0.00001]; %Desplazamiento, Velocidad, Angulo, Velocidad angular, Integrador
QQ = diag(dd);
RR = 1e-5;                    

KK = dlqr(AA,BB,QQ,RR);
K = KK(1:4);
KI = -KK(5);

% Controlador  de 10 a 0 m
%matrices para m*10
m2 = m*10;

Ac_m2=[0       1               0       0;  %x1=delta - desplazamiento
    0    -Fricc/M         -m2*g/M    0;  %x2=delta_p
    0       0               0       1;  %x3=phi - angulo
    0  -Fricc/(l*M)  -g*(m2+M)/(l*M) 0];  %x4=phi_p

sys_c_m2 = ss(Ac_m2, Bc, Cc, Dc);
sys_d_m2=c2d(sys_c_m2,Ts,'zoh'); 

%Matrices de lazo abierto en el equilibrio estable
A_m2 = sys_d_m2.a;
B_m2 = sys_d_m2.b;
C_m2 = sys_d_m2.c;

AA_m2=[A_m2,zeros(4,1);-Cref*A_m2,eye(1)];% para el integrador de m2

% parametros DLQR
dd_m2 = [.1 1e-3 1e-3 .1 0.001]; %Desplazamiento, Velocidad, Angulo, Velocidad angular, Integrador
QQ_m2 = diag(dd);
RR_m2 = .008;                    

KK_m2 = dlqr(AA_m2,BB,QQ_m2,RR_m2);
K_m2 = KK_m2(1:4);
KI_m2 = -KK_m2(5);




%Observador
Ao = A';
Bo = Cc';
Co = B';

%parametros DLRQ- Observador
do = [0.001 1000 0.5 0.0001]; %Desplazamiento, Velocidad, Ã?ngulo, Velocidad angular
Qo = diag(do); 
Ro = diag([80 10000]);

Kko = dlqr(Ao,Bo,Qo,Ro);
Ko=Kko';
t = 0;

%Ganancia de prealimentacion
Gj = inv(Cref*inv(eye(4)-A+B*K)*B);

x = [0;0;alfa(1);0];
p = x(1);
p_p = x(2);
alfa = x(3);
w = x(4);

tita_pp(1) = 0;
h = Ts/20;
u = []; 
i = 1;
u_k(1) = 0; 
ref = 10; 
flag = 0;
v(1) = 0;
x = [0; 0; alfa(1);0];
x_hat = [0;0;pi;0];
xop=[0 0 pi 0]';

for ki=1:Kmax
        
    y_sal=Cc*x;  %Salida de dos componentes
    y_sal_o=Cc*(x_hat-xop);
    v(ki+1)=v(ki)+ref-y_sal(1);
    
    %Ley de control
    u1(ki)=-K*(x-xop)+KI*v(ki+1); %color = '';%Sin observador
    %u1(ki)=-K*x_hat+KI*v(ki+1);color = ''; %Con observador
    
    %Zona Muerta
    zona_muerta=0;
    if(abs(u1(ki))<zona_muerta)
        u1(ki)=0;               
    else
        u1(ki)=sign(u1(ki))*(abs(u1(ki))-zona_muerta);
    end
    %-----------------------------------------------------
    
    for kii=1:Ts/h
        
        u(i)=u1(ki);
        p_pp=(1/(M+m))*(u(i)-m*l*tita_pp*cos(alfa(i))+m*l*w(i)^2*sin(alfa(i))-Fricc*p_p(i));
        tita_pp=(1/l)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
        p_p(i+1)=p_p(i)+h*p_pp;
        p(i+1)=p(i)+h*p_p(i);
        w(i+1)=w(i)+h*tita_pp;
        alfa(i+1)=alfa(i)+h*w(i);
        if(p(i)>=9.99)
            if(flag==0)
                ref=0
                m=m*10;
                flag=1;
                K=K_m2;
                KI=KI_m2;
            end
        end
        i=i+1;
    end
    x=[p(i-1); p_p(i-1); alfa(i-1); w(i-1)];
    x_hat=A*x_hat+B*u1(ki)+Ko*(y_sal-y_sal_o)+xop;
end

u(i)=u1(ki);
t=0:h:T;

figure(1);
subplot(3,2,1); grid on; hold on;
plot(t,w,color,'LineWidth',1.5);grid on; title('Velocidad angular \omega');

subplot(3,2,2); grid on; hold on;
plot(t,alfa/pi*100,color,'LineWidth',1.5); title('Ángulo \phi');xlabel('Tiempo');

subplot(3,2,3); grid on; hold on;
plot(t,p,color,'LineWidth',1.5);title('Posición grúa \theta');xlabel('Tiempo');

subplot(3,2,4); grid on; hold on;
plot(t,p_p,color,'LineWidth',1.5);title('Velocidad de grúa \theta_p');

subplot(3,1,3); grid on; hold on;
plot(t,u,color,'LineWidth',1.5);title('Acción de control u');xlabel('Tiempo en Seg.');
% 
% figure(2);
% subplot(2,1,1);grid on; hold on;
% plot(alfa,w,'');
% title('Ángulo vs Velocidad angular');
% xlabel('Ángulo');ylabel('Velocidad angular');
% 
% subplot(2,1,2);grid on; hold on;
% plot(p,p_p, '');
% title('Distancia vs velocidad');
% xlabel('Distancia');ylabel('Velocidad');

