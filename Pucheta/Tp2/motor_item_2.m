% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 2 - Caso de estudio 1 - 
%   Inciso 2  
% Considerar que no puede medirse la corriente y sólo pueda medirse el ángulo, 
% por lo que debe implementarse un observador. Obtener la simulación en las mismas 
% condiciones que en el punto anterior, y superponer las gráficas para comparar.
% 
%%

clc;clear all;close all;
%Parametros de motor
Laa=5e-3;J=0.004;
Ra=0.2;Bm=0.005;
Ki=6.5e-5;
Km=0.055;

%Matrices de estado
A=[-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B=[1/Laa; 0; 0];
C=[0 0 1];%salida solo el angulo
D=[0];


%Matrices ampliadas
Aa=[A zeros(3,1); -C 0];
Ba=[B(:,1); 0];
Ca=[C 0];

%LQR
Q=diag([10 1/100 1/100 10000/2.5]);    R=10;
Kamp=lqr(Aa,Ba,Q,R);

%Observador--------------------------------------
Ao=A';
Bo=C';
Co=B';

% Qo=diag([.0001 .1 .01]);    Ro=1000;
Qo=1e2*diag([40000 0 0]);    Ro=.01;
Ko=lqr(Ao,Bo,Qo,Ro);

Tf=100;h=1e-4;pasos=Tf/h; t=linspace(0,Tf,pasos); %Con Tf=300 y tc=100
Ref=linspace(0,0,pasos);                    % se logra una mejora
TL=linspace(0,0,pasos);                     % en la respuesta
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

estados=[ia(1); omega(1); theta(1)];
estados_obs=[ia(1);omega(1);theta(1)];

x=[ia(1) omega(1) theta(1)]';
xobs=[0 0 0]'; %inicializacion para el observador

psi(1)=0;
integracion(1)=psi(1);

for i=1:round(Tf/h)
    
    psi_p=Ref(i)-Ca(1:3)*estados-Ca(4)*integracion;
    psi(i)=integracion+psi_p*h;
    
    u(i) = -Kamp(1:3)*estados_obs-Kamp(4)*psi(i);
    %Variables del sistema lineal
    ia(i)= x(1);
    omega(i)= x(2);
    theta(i)= x(3);
    
    x1_p=-Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
    x2_p=Ki*x(1)/J-Bm*x(2)/J-TL(i)/J;
    x3_p=x(2);

    xp=[x1_p; x2_p; x3_p];
    x=x+h*xp;
    
    %------con Observador----------------------
  
    ia_0(i)= xobs(1);
    omega_0(i)= xobs(2);
    theta_0(i)= xobs(3);
    
    y_sal_o(i) = C * estados_obs;
    y_sal(i)   = Ca(1:3) * estados + Ca(4)*integracion;
    x_hat_p     = A*xobs+B*u(i)+Ko*(y_sal(:,i)-y_sal_o(:,i));
    xobs       = xobs + x_hat_p*h;
    %--------------------------------------------
    
    estados=[ia(i);omega(i);theta(i)];
    integracion=psi(i);
    estados_obs=[ia_0(i);omega_0(i);theta_0(i)];
    
end

%------Sistema sin observador---------------------------------------------

%condiciones iniciales
ia_so(1)=0;                %Corriente de armadura x1
theta_so(1)=0;             %Valocidad angular x2
omega_so(1)=0;             %Posicion angular x3

estados_so=[ia_so(1);omega_so(1); theta_so(1)];
x_so=[ia_so(1) omega_so(1) theta_so(1)]';

psi_so(1)=0;
integracion_so(1)=psi_so(1);

for i=1:round(Tf/h)
    
    psi_p_so=Ref(i)-C*estados_so;
    psi_so(i)=integracion_so+psi_p_so*h;
    
    u_so(i) = -Kamp(1:3)*estados_so-Kamp(4)*psi_so(i);
    %Variables del sistema lineal
    ia_so(i)= x_so(1);
    omega_so(i)= x_so(2);
    theta_so(i)= x_so(3);
    
    x1_p=-Ra*x_so(1)/Laa-Km*x_so(2)/Laa+u_so(i)/Laa;
    x2_p=Ki*x_so(1)/J-Bm*x_so(2)/J-TL(i)/J;
    x3_p=x_so(2);
    
    xp_so=[x1_p; x2_p; x3_p];
    x_so=x_so+h*xp_so;
    
    estados_so=[ia_so(i);omega_so(i);theta_so(i)];
    integracion_so=psi_so(i);
end

%-------------------------------------------------------------------------


subplot(2, 2, 1);
hold on
plot(t,ia,'b');
plot(t,ia_so,'r');
hold off
legend({'Con observador','Sin observador'},'Location','southeast')
title('Corriente de armadura i_a');
xlabel('Tiempo (seg.)');
ylabel('Corriente (A)');
grid on;

subplot(2, 2, 2);
hold on
plot(t,omega,'b');
plot(t,omega_so,'r');
hold off
legend({'Con observador','Sin observador'})
title('Velocidad angular \omega_r');
xlabel('Tiempo (seg.)');
ylabel('Velocidad angular (rad/s)');
grid on;

subplot(2, 2, 3);
hold on
plot(t,theta,'b');
plot(t,theta_so,'r');
plot(t,Ref,'m');
hold off
legend({'Con observador','Sin observador','Referencia'})
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(2, 2, 4);
hold on
plot(t,u,'b');
plot(t,u_so,'r');
hold off
legend({'Con observador','Sin observador'})
title('Accion de control u_t');
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;