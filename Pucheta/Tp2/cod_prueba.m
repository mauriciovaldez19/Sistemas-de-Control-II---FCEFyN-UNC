clc; clear;
%Parametros y variables de simulacion
%LAA=366e-6; J=5e-9; RA=55.6; B=0; Ki=6.49e-3; Km=6.53e-3;
LAA=5e-3; J=0.004; RA=0.2; B=0.005; Ki=6.5e-5;Km=0.055;
Tf=12;h=1e-4; pasos=Tf/h; t=0:h:Tf;
ia=0:h:Tf; w_p=0:h:Tf; w=0:h:Tf; theta=0:h:Tf; t=0:h:Tf;
u =linspace(0,Tf,pasos+1);
ref=linspace(0,Tf,pasos+1);
%Definicion variables
tRef=(pi/2);
Tl=1.15e-7;
tc=2;
est=0;
ii=1;
kk=0;
%Condiciones iniciales
ia(1)=0; w(1)=0; theta(1)=0; w_p(1)=0; u(1)=0; TL=Tl;
%Matrices de estados
Mat_A=[ -RA/LAA -Km/LAA 0; %x1=corriente
 Ki/J -B/J 0; %x2=velocidad angular
 0 1 0]; %x3=posicion
Mat_B=[ 1/LAA 0 ; 
 0 -1/J;
 0 0];
Mat_C=[0 0 1]; %La salida monovariable es posición y ángulo
Mat_M=[Mat_B(:,1) Mat_A*Mat_B(:,1) Mat_A^2*Mat_B(:,1)];%Matriz Controlabilidad
%Cálculo del controlador por asignación de polos
auto_val=eig(Mat_A);
pol_caract_A=poly(Mat_A);
Mat_W=[ pol_caract_A(3) pol_caract_A(2) 1;
 pol_caract_A(2) 1 0;
 1 0 0];
 
Mat_T=Mat_M*Mat_W;
A_controlable=inv(Mat_T)*Mat_A*Mat_T; %Verificación de que T esté bien
 %Si la inversa existe es controlable
%Ubicación de los polos de lazo cerrado
P1=-1e5;P2=-80; P3=-50;
alfa_i=poly([P1 P2 P3]); %Genero la ecuacion de los polos deseada 
K=(fliplr(alfa_i(2:4)- pol_caract_A(2:4)))*inv((Mat_T));
eig(Mat_A-Mat_B(:,1)*K); %Verifico que todos los polos estén en el semiplano izquierdo
G=-inv(Mat_C*inv(Mat_A-(Mat_B(:,1)*K))*Mat_B(:,1)); %Matriz G para correr la referencia del 0
%Observador
Mat_AO=Mat_A';
Mat_BO=Mat_B(:,1)'; 
Mat_CO=Mat_C';
Mat_M_Dual=[Mat_CO Mat_AO*Mat_CO Mat_AO^2*Mat_CO];%Matriz Controlabilidad
 % del observador
% Ubicacion del Observador
% Algunas veces más rápido que el controlador
alfa_io=poly([P1 P2 P3]*150);
Mat_TO=Mat_M_Dual*Mat_W;
KO=(fliplr(alfa_io(2:end)-pol_caract_A(2:end))*inv(Mat_TO))';
eig(Mat_AO'-KO*Mat_C); %Verifico que todos los polos estén en el semiplano izquierdo
x_hat=[0 0 0]'; %Inicializo el Observador
%Simulacion
while(ii<(pasos+1))
 kk=kk+h;
 if(kk>tc)
 tRef=tRef*(-1);
 if(est==0)
 TL=0; %Establece la referencia a seguir
 est=1;
 else
 TL=Tl;
 est=0;
 end
 kk=0;
 end
 
 ref(ii)=tRef;
 estado=[ia(ii); w(ii); theta(ii)];
 
 %Ley de control
 u(ii)= -K*estado+G*ref(ii); color='b'; %Sin observador
% u(ii)= -K*x_hat(:,1)+G*ref(ii); color='r'; %Con observador
 U=u(ii);
 %Satuador
 if U>12
 u(ii)=12;
 end
 
 if U<-12
 u(ii)=-12;
 end 
 %Integracion Euler
 w_pp =(-w_p(ii)*(RA*J+LAA*B)-w(ii)*(RA*B+Ki*Km)+u(ii)*Ki)/(J*LAA);
 ia_p=(1/LAA)*(-RA*ia(ii)-Km*w(ii)+u(ii));
 w_p(ii+1)=w_p(ii)+h*w_pp-(1/J)*TL;
 ia(ii+1)=ia(ii)+h*ia_p;
 w(ii+1)=w(ii)+h*w_p(ii);
 theta(ii+1)=theta(ii)+h*w(ii);
 
 
 y_salO(ii)=Mat_C*x_hat(:,1);
 y_sal(ii)=Mat_C*estado;
 x_hatp=Mat_A*x_hat+Mat_B*u(ii)+KO*(y_sal(ii)-y_salO(ii));
 x_hat=x_hat+h*x_hatp;
 
 ii=ii+1;
 
end
%Graficos
figure(1);hold on ;
subplot(2,2,1);plot(t,ia,color);grid on; title('i_a , Corriente [A]');hold on; xlabel('Tiempo en Seg.');
subplot(2,2,2);plot(t,theta,color);grid on;hold on;
subplot(2,2,2);plot(t,ref,'k');grid on;title('\theta_t , Posicion angular [rad]');hold on; xlabel('Tiempo en Seg.');
subplot(2,1,2);plot(t,u,color);grid on;title('u_t , Acción de control [V]');xlabel('Tiempo en Seg.');
