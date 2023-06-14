%Caso de estudio 2. Sistema lineal de cuatro variables de estado2
clear; clc;

%Se declaran los parámetros del sistema
w = 9;          %Frecuencia natural
a = 0.07;
b = 5;
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
       
Cc=[0 0 0 1]   %La salida es la altura

Dc=[0]

%val=real(1/eig(Ac))

%Punto N°3: Establecer el Ts más adecuado
%Tiempos
Ts = 0.001;             %Tiempo de muestreo  
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

%Pasaje de tiempo continuo a tiempo discreto
sys_c = ss(Ac, Bc, Cc, Dc);
sys_d=c2d(sys_c,Ts,'zoh');

%TIEMPO DISCRETO
%Matrices de lazo abierto
A = sys_d.a;
B = sys_d.b;
C = sys_d.c;

%Punto N°1: se implementa un controlador en  tiempo discreto

%DISEÑO DEL CONTROLADOR
%Método: DLQR (Discret Linear Quadratic Regulator)

%Se verifica la controlabilidad
%Es alcanzable <-> rank(Ma)= 4, ya que n=4
%Es controlable <-> rank(Ma)= rank(Mc)
 
%Matriz de Ma y Mc
Ma = [B A*B A^2*B A^3*B];       %Alcanzabilidad
Mc = [B A*B A^2*B A^3*B A^4];   %Controlabilidad
 
Alcanzabilidad = rank(Ma)
Controlabilidad = rank(Mc)

%Parámetros del controlador DLQR
d = [1 1000 150000 0.05]; %Alpha, Phi, Phi_p, Altura
Q = diag(d);
R = 1;                    %Dimensiona la acción de control

K = dlqr(A, B, Q, R)

%Punto N°2: se implementa un controlador con observador
%Observador
Ao = A'
Bo = Cc'
Co = B'

%Parámetros del controlador DLQR observador
do = [0.0000001 0.0000001 0.000001 100];
Qo = diag(do); 
Ro = 100;

Ko = (dlqr(Ao,Bo,Qo,Ro))'

%Ganancia de prealimentación
Gj = inv(Cc*inv(eye(4)-A+B*K)*B);

%Se definen los integradores y se setean a cero
alfap = 0:At:T;
fip = 0:At:T;
fipp = 0:At:T;
hp = 0:At:T;

x_hat = [0; 0 ;0; 0];
i=1;

%Simulación
while(i<(pasos+1))
    x=[alfa(i);fi(i);fip(i);h(i)];
    
    %Ley de control
    u(i)=-K*x+Gj*href; %Sin observador
    %u(i) = -K*x_hat+Gj*href; %Con observador
    
    %Punto N°4: se aumenta la zona muerta hasta alcanzar valores inaceptables
    zona_muerta=0;
    
    %Zona Muerta
    if(abs(u(i))<zona_muerta)
        u(i)=0;
    else
        u(i)=sign(u(i))*(abs(u(i))-zona_muerta);
    end
    
    %-------------------------------------------
    
    %EDOS, Sistema Real:
    alfap = a*(fi(i)-alfa(i));
    fipp = -w^2*(fi(i) - alfa(i) - b*u(i));
    hp = c*alfa(i);
    y_sal = Cc*x;
    y_sal_o=Cc*x_hat;
    x_hat = A*x_hat+B*u(i)+Ko*(y_sal - y_sal_o);
    
    %Integración por Euler
    alfa(i+1) = alfa(i) + alfap*At;
    fip(i+1)=fip(i)+fipp*At;
    fi(i+1) = fi(i) + fip(i)*At;
    h(i+1) = h(i) + hp*At;
    y_sal(i)=Cc*x;
    i=i+1;
    
end

color='';

%Gráficos
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
plot(t,u,color);grid on;title('Acción de control u');xlabel('Tiempo');
