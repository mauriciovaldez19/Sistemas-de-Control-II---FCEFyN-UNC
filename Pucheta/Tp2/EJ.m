% Función que representa la ecuación diferencial
function dy = ecuacion_diferencial(t, y)
    % Definir aquí la ecuación diferencial a resolver
    % Por ejemplo: dy = -2 * t * y;
    dy = -2 * t * y;
end

% Parámetros de la integración
t0 = 0;         % Valor inicial de t
y0 = 1;         % Valor inicial de y
h = 0.1;        % Tamaño del paso de integración
tf = 1;         % Valor final de t

% Cálculo de la solución mediante el método de Euler
t = t0:h:tf;    % Vector de valores de t
n = length(t);  % Número de pasos de integración
y = zeros(size(t));   % Vector de valores de y
y(1) = y0;      % Valor inicial de y

for i = 1:(n-1)
    dy = ecuacion_diferencial(t(i), y(i));    % Calcular la derivada en el punto actual
    y(i+1) = y(i) + h * dy;  % Actualizar el valor de y mediante el método de Euler
end

% Graficar la solución
plot(t, y, 'b-o');
xlabel('t');
ylabel('y');
title('Método de integración por Euler');
grid on;
