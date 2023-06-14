%Caso de estudio 2. Sistema lineal de cuatro variables de estado2
clear; clc;

%Se declaran los parámetros del sistema
w = 9;          %Frecuencia natural
a = 0.07;       %cte positiva
b = 5;          %cte positiva
c = 150;        %Velocidad m/s

%TIEMPO CONTINUO
%Matrices de lazo abierto
Ac=[-a    a   0  0;          %x1=Alpha
     0    0   1  0;          %x2=phi
    w^2 -w^2  0  0;          %x3=phi_p
     c    0   0  0]          %x4=h
       
Bc=[   0;
       0;
    (w^2)*b;
       0]
       
Cc=[0 0 0 1;0 1 0 0]   %La salida es la altura y angulo

Dc=[0]

%val=real(1/eig(Ac))

% Polos deseados para el controlador
p1 = -15+15i;
p2 = -15-15i;
p3 = -0.5+0.5i;
p4 = -0.5-0.5i;
polos = [p1 p2 p3 p4];

Kk=place(Ac,Bc,polos)


% %LQR
% d = [.1 .1 1 1]; %Alpha, Phi, Phi_p, Altura
% Q = diag(d);
% R = 1000;                  
% Kk = lqr(Ac, Bc, Q, R)

%Punto N°2: se implementa un controlador con observador
%Observador
Ao = Ac'
Bo = Cc'
Co = Bc'


do = [1 1 1 1];
Qo = diag(do); 
Ro=[1e2    0  ;
      0    1e2];;

Ko = (lqr(Ao,Bo,Qo,Ro))'

%Ganancia de prealimentación
Gj = inv(Cc*inv(eye(4)-Ac+Bc*Kk)*Bc);
T = 70;                 %Tiempo de simulación
At = 0.001;             %Tiempo de integración 
pasos = T/At; 

%Asignaciones
t = 0:At:T;
alfa = zeros(1,pasos);
fi = zeros(1,pasos);
fip = zeros(1,pasos);
h = zeros(1,pasos); 
u = linspace(0,0,pasos+1);

%Condiciones Iniciales
alfa(1) = 0; 
fi(1) = 0;
fip(1) = 0;
h(1) = -500;            %Altura inicial - Modificar a 500
u(1) = 0;
href = 100;             %Altura de referencia - Modificar a -100


%Se definen los integradores y se setean a cero
alfap = 0:At:T;
fip = 0:At:T;
fipp = 0:At:T;
hp = 0:At:T;

x_hat = [0; 0 ;0; 0];
i=1;

%Simulaci�n
while(i<(pasos+1))
    x=[alfa(i);fi(i);fip(i);h(i)];
    
    %Ley de control
    u(i)=-Kk*x+Gj*href; %Sin observador
    %u(i) = -K*x_hat+Gj*href; %Con observador
    
    
    %EDOS, Sistema Real:
    alfap = a*(fi(i)-alfa(i));
    fipp = -w^2*(fi(i) - alfa(i) - b*u(i));
    hp = c*alfa(i);
    y_sal = Cc*x;
    y_sal_o=Cc*x_hat;
    x_hat = Ac*x_hat+Bc*u(i)+Ko*(y_sal - y_sal_o);
    
    %Integraci�n por Euler
    alfa(i+1) = alfa(i) + alfap*At;
    fip(i+1)=fip(i)+fipp*At;
    fi(i+1) = fi(i) + fip(i)*At;
    h(i+1) = h(i) + hp*At;
    y_sal(i)=Cc*x;
    i=i+1;
    
end

color='';

%Gr�ficos
figure(1);hold on;
subplot(3,2,1);hold on;grid on;
plot(t,h,color);hold on;title('Altura h');xlabel('Tiempo');

subplot(3,2,2);hold on;grid on;
plot(t,alfa,color);title('\alpha');xlabel('Tiempo');

subplot(3,2,3);hold on;grid on;
plot(t,fi,color);title('\phi');xlabel('Tiempo');

subplot(3,2,4);hold on; grid on;
plot(t,fip, color);title('\phi_p');xlabel('Tiempo');

subplot(3,1,3);hold on;grid on;
plot(t,u,color);grid on;title('Acci�n de control u');xlabel('Tiempo');
