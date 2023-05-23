%%Tp2 - Ing Laboret
%Alumno: Valdez Benavidez, Mauricio Luciano

% Datos: Polo1=-2       Polo2=0     Cero=-10
%Ganancia=10    Sobrepaso=5    tiempo2%=3   error=0 periodo=0.15
%% Definicion FdT continua
G=zpk([-10],[-2 0],[10]);

%%  Definicion FdT discreta
T=0.15;
Gd=c2d(G,T,'zoh');
%% Obtención de W0;Wd;td;psita
S=5;    tr=3;   error=0;    T=0.15;

psita=(-log(S/100))/sqrt(pi^2+log(S/100)^2)
W0=4/(tr*psita)
Wd=W0*sqrt(1-psita^2)
td=(2*pi)/Wd
%% Obtención muestras por ciclo de la frec. amortiguada wd
m=td/T
%% Ubicacion de los polos en el plano z
r=exp(-psita*W0*T)
ang=rad2deg(Wd*T)

rect=r*(cos(ang)+j*sin(ang))

%% diseño de controlador con sisotool
%  sisotool(Gd)
% 
% %%Verificacion polos, ceros y respuesta al escalon temporal
% C %muestra el compensador importado de sisotool
% F=feedback(C*Gd,1) % sistema de lazo cerrado
% pole(F)
% zero(F)
% pzmap(F)
% step(F) % respuesta al escalon
%pid(C) %te da los valores de los K individuales