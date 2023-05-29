% %%
% SE MODIFICO EL CODIGO DEL EJEMPLO 5.2 PARA LA CONSIGNA DEL TP2
% CASO 2 ITEM 3 Y POSIBLE 4
%%
clc;clear all;
m=.1;Fricc=0.1; long=1.6;g=9.8;M=1.5;
h=1e-3;tiempo=(20/h);p_pp=0;tita_pp=0; t=0:h:tiempo*h;
omega=0:h:tiempo*h; alfa=0:h:tiempo*h; p=0:h:tiempo*h;
p_p=0:h:tiempo*h; u=linspace(0,0,tiempo+1);
%Condiciones iniciales
alfa(1)=.1; color='r';
% alfa(1)=.2; color='g';
% alfa(1)=.7; color='b';
ref=-10;
omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;indice=0;
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0]
Mat_B=[0; 1/M; 0; -1/(long*M)]
Mat_C=[1 0 1 0]; %La salida es posición y ángulo
% % % Construcción del sistema ampliado
% % Mat_Aa=[Mat_A zeros(4,1);-Mat_C 0];
% % Mat_Ba=[Mat_B;0];
% % Mat_Ma=[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba Mat_Aa^4*Mat_Ba];%Matriz Controlabilidad
% % %Cálculo del controlador por asignación de polos
% % auto_val=eig(Mat_Aa);
% % c_ai=conv(conv(conv(conv([1 -auto_val(1)],[1 -auto_val(2)]),[1 -auto_val(3)]),[1 -auto_val(4)]),[1 -auto_val(5)]);
% % Mat_Wa=[c_ai(5) c_ai(4) c_ai(3) c_ai(2) 1;c_ai(4) c_ai(3) c_ai(2) 1 0;c_ai(3) c_ai(2) 1 0 0;c_ai(2) 1 0 0 0;1 0 0 0 0];
% % Mat_Ta=Mat_Ma*Mat_Wa;
% % A_controlable=inv(Mat_Ta)*Mat_Aa*Mat_Ta; %Verificación de que T esté bien
% % %Ubicación de los polos de lazo cerrado en mui:
% % mui(1)=-.7;mui(2)=-.7; mui(3)=-10 + 0.4i;mui(4)=conj(mui(3));mui(5)=-1;
% % alfa_ia=conv(conv(conv(conv([1 -mui(3)],[1 -mui(4)]),[1 -mui(2)]),[1 -mui(1)]),[1 -mui(5)]);
% % Ka=(alfa_ia(2:6)-c_ai(2:6))*inv(Mat_Ta);
% % eig(Mat_Aa-Mat_Ba*Ka)
% % K=Ka(1:4); KI=-Ka(5); %Los valores del controlador de obtienen del K ampliado
% % psi(1)=0;
%METODO LQR
D = [1e2 100 1e4 20];  %Velocidad angular, Posicion angular, Posicion Carrito, Velocidad Carrito 
Q = diag(D);
R = 1000;
K = lqr(Mat_A,Mat_B,Q,R);
KI=-inv(Mat_C*(eye(4)-inv(Mat_A+Mat_B*K))*Mat_B); %%Es igual a Klqr(1)



while(i<(tiempo+1))
 estado=[p(i); p_p(i); alfa(i); omega(i)];
%  psi_p=ref-Mat_C* estado;
%  psi(i+1)=psi(i)+psi_p*h;
  u(i)=-K*estado+K(1)*ref; %+KI*psi(i+1);
 p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
 tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
 p_p(i+1)=p_p(i)+h*p_pp;
 p(i+1)=p(i)+h*p_p(i);
 omega(i+1)=omega(i)+h*tita_pp;
 alfa(i+1)=alfa(i)+h*omega(i);
 y_sal(i)=Mat_C*estado; i=i+1;
end
figure(1);hold on;
subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad ángulo');hold on;
subplot(3,2,2);plot(t,alfa,color);grid on;title('Ángulo');hold on;
subplot(3,2,3); plot(t,p,color);grid on;title('Posición carro');hold on;
subplot(3,2,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;
figure(2);hold on;subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;