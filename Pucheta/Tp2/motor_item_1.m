clc; clear ; close all
%EN ESTE CASO VAMOS A AGREGAR UN INTEGRADOR Y UNA GANANCIA K_I
%DEFINO PARAMETROS

LAA=5e-3; 
J=0.004; 
RA=0.2; 
Bm=0.005; 
Ki=6.5e-5;
Km=0.055;

%DEFINO MATRICES
%X=[ia ; tita ; w];
A=[-RA/LAA 0 -Km/LAA  ; 0 0 1 ; Ki/J 0 -Bm/J];
B=[1/LAA; 0; 0];
C=[0 1 0];
D=[0];
%pero ahora como vamos a agregar integrador las matrices A y B se ven
%modificadas
An=[A zeros(3,1); -C 0];
Bn=[B ; 0];
Cn=[1 0];

%implementacion de funciones a usar
tf=100; dt=1*10^-4; t=0:dt:(tf-dt); periodo=40;%[seg]
torq=1.15*10^-3;

Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3

%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
Q=diag([1/10 1/10 1 20]); R=1/10;
Ka=lqr(An,Bn,Q,R);
K_i= -Ka(4);
K=Ka(1:3);

Kn=[K -K_i];

%iteracion
n=round(tf/dt);
X=zeros(3,n);
X(1,1)=0; %ia inicial
X(2,1)=0; %tita inicial
X(3,1)=0; %w inicial
psi(1)=0; %psi inicial
%DEFINO CONDICIONES INICIALES

for i=1:1:n-1
    X_a=[X(1,i); X(2,i) ; X(3,i)];%[ia ; tita ; w]
    psi_p=Ref(i)-C*X_a;
    psi(i+1)=psi(i)+psi_p*dt;
    U=-K*X_a+K_i*psi(i+1);
 if U>12
 U=12;
 end
 
 if U<-12
 U=-12;
 end 
    
    Xp_1=-RA/LAA*X_a(1)-Km/LAA*X_a(3)+1/LAA*U;  %ia_p
    Xp_2= X_a(3);                               %tita_p
    Xp_3=Ki/J*X_a(1)-Bm/J*X_a(3)-1/J*TL(i);%    %W_p
    
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3];
    
    Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    
    X(1,i+1)=Xf(1);
    X(2,i+1)=Xf(2);
    X(3,i+1)=Xf(3);
    
end
%ploteo de entrada con ganancia de prealimentacion U  y perturbacion TL
% figure(2)
% subplot(2,1,1)
% hold on; grid on;
% plot(Ref);title('Referencia de Entrada');xlabel('tiempo[s]');ylabel('angulo');
% subplot(2,1,2)
% hold on;grid on;
% plot(TL);title('Torque de perturbación');xlabel('Tiempo');ylabel('Torque');


figure(1)
subplot(2,1,1);
plot(t,Ref);
grid on
hold on
plot(t,X(2,:),'r');title('angulo tita con integrador y sin observador ');xlabel('tiempo[s]');ylabel('angulo[rad]');
legend('REF','var de estado')
subplot(2,1,2);
%plot(t,Ref);
grid on
hold on
plot(t,X(1,:),'r');title('corriente i_a con integrador y sin observador ');xlabel('tiempo[s]');ylabel('corriente[A]');
legend('var de estado')




