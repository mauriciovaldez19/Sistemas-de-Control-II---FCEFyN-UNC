% Funci�n que representa la ecuaci�n diferencial
function dy = ecuacion_diferencial(t, y)
    % Definir aqu� la ecuaci�n diferencial a resolver
    % Por ejemplo: dy = -2 * t * y;
    dy = -2 * t * y;
end

% Par�metros de la integraci�n
t0 = 0;         % Valor inicial de t
y0 = 1;         % Valor inicial de y
h = 0.1;        % Tama�o del paso de integraci�n
tf = 1;         % Valor final de t

% C�lculo de la soluci�n mediante el m�todo de Euler
t = t0:h:tf;    % Vector de valores de t
n = length(t);  % N�mero de pasos de integraci�n
y = zeros(size(t));   % Vector de valores de y
y(1) = y0;      % Valor inicial de y

for i = 1:(n-1)
    dy = ecuacion_diferencial(t(i), y(i));    % Calcular la derivada en el punto actual
    y(i+1) = y(i) + h * dy;  % Actualizar el valor de y mediante el m�todo de Euler
end

% Graficar la soluci�n
plot(t, y, 'b-o');
xlabel('t');
ylabel('y');
title('M�todo de integraci�n por Euler');
grid on;
