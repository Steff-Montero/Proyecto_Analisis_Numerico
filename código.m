% Parámetros del sistema
m = 1.0;                           % Masa(kg)
c = 0.5;                            % Coeficiente de amortiguamiento
k = 4.0;                            % Constante del resorte
F0 = 1.0;                          % Amplitud de la fuerza externa
omega = 2.0;                    % Frecuencia de la fuerza externa

% Definimos el sistema de ecuaciones diferenciales de primer orden
sistema = @(t, estado) [estado(2); (F0 * cos(omega * t) - c * estado(2) - k * estado(1)) / m];

% Condiciones iniciales
x0 = 0.0;                           % Posición inicial
v0 = 0.0;                           % Velocidad inicial
estado_inicial = [x0; v0];    
t0 = 0.0;                           % Tiempo inicial
t_max = 20.0;                   % Tiempo final
h = 0.01;                       % Tamaño del paso

% Inicializamos vectores para almacenar los resultados
t = t0:h:t_max;                 % Tiempos
N = length(t);
estado = zeros(2, N);           % Cada fila representa x y v
estado(:, 1) = estado_inicial;

% Método de Runge-Kutta de cuarto orden
for i = 1:N-1
    k1 = sistema(t(i), estado(:, i));
    k2 = sistema(t(i) + h/2, estado(:, i) + h * k1 / 2);
    k3 = sistema(t(i) + h/2, estado(:, i) + h * k2 / 2);
    k4 = sistema(t(i) + h, estado(:, i) + h * k3);
    estado(:, i+1) = estado(:, i) + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end

% Extraemos el desplazamiento x(t) y la velocidad v(t)
x = estado(1, :);
v = estado(2, :);

% Graficamos el desplazamiento y la velocidad
figure;

% Gráfica del desplazamiento
subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 1.5);
grid on;
title('Desplazamiento $x(t)$', 'Interpreter', 'Latex');
xlabel('Tiempo (s)');
ylabel('Desplazamiento (m)');
legend('x(t)', 'Interpreter', 'Latex');

% Gráfica de la velocidad
subplot(2, 1, 2);
plot(t, v, 'r', 'LineWidth', 1.5);
grid on;
title('Velocidad $v(t)$', 'Interpreter', 'Latex');
xlabel('Tiempo (s)');
ylabel('Velocidad (m/s)');
legend('v(t)', 'Interpreter', 'Latex');

%% Vibración amortiguada
% Parámetros del sistema
m = 1.5;      % Masa (kg)
k = 9.0;      % Constante del resorte (N/m)

% Calcular amortiguamiento crítico
c_critico = 2 * sqrt(m * k);

% Valores de amortiguamiento para cada régimen
c_values = [0.1 * c_critico, c_critico, 4 * c_critico]; % Subamortiguado, crítico, sobreamortiguado
labels = {'Subamortiguado', 'Crítico', 'Sobreamortiguado'};

% Condiciones iniciales
x0 = 0.2; % Posición inicial (m)
v0 = -0.5; % Velocidad inicial (m/s)
estado_inicial = [x0; v0];
t0 = 0.0; % Tiempo inicial
t_max = 10.0; % Tiempo final
h = 0.01; % Tamaño del paso

% Inicializar gráficos para desplazamiento y energía
figure;
subplot(2, 1, 1); % Subgráfico para desplazamiento
hold on;
title('Regímenes de amortiguamiento: Desplazamiento');
xlabel('Tiempo (s)');
ylabel('Desplazamiento (m)');
grid on;

subplot(2, 1, 2); % Subgráfico para energía total
hold on;
title('Efecto del amortiguamiento en la energía total');
xlabel('Tiempo (s)');
ylabel('Energía total (J)');
grid on;

% Resolver para cada valor de c
for idx = 1:length(c_values)
    c = c_values(idx);
    
    % Definir el sistema de ecuaciones diferenciales
    sistema = @(t, estado) [ ...
        estado(2); % dx/dt = v
        (-c * estado(2) - k * estado(1)) / m % dv/dt
    ];
    
    % Inicializar vectores para almacenar resultados
    t = t0:h:t_max;
    N = length(t);
    estado = zeros(2, N);
    estado(:, 1) = estado_inicial;
    
    % Método de Runge-Kutta de cuarto orden (RK4)
    for i = 1:N-1
        k1 = sistema(t(i), estado(:, i));
        k2 = sistema(t(i) + h/2, estado(:, i) + h * k1 / 2);
        k3 = sistema(t(i) + h/2, estado(:, i) + h * k2 / 2);
        k4 = sistema(t(i) + h, estado(:, i) + h * k3);
        estado(:, i+1) = estado(:, i) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
    
    % Extraer el desplazamiento y velocidad
    x = estado(1, :);
    v = estado(2, :);
    
    % Calcular energía cinética, potencial y total
    energia_cinetica = 0.5 * m * v.^2;
    energia_potencial = 0.5 * k * x.^2;
    energia_total = energia_cinetica + energia_potencial;
    
    % Graficar desplazamiento
    subplot(2, 1, 1);
    plot(t, x, 'LineWidth', 1.5, 'DisplayName', labels{idx});
    
    % Graficar energía total
    subplot(2, 1, 2);
    plot(t, energia_total, 'LineWidth', 1.5, 'DisplayName', labels{idx});
end

% Mostrar leyendas
subplot(2, 1, 1);
legend show;
subplot(2, 1, 2);
legend show;
hold off;
