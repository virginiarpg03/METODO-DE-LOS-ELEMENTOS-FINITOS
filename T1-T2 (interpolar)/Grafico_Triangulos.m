function [] = Grafico_Triangulos(f)

    [xi, yi] = meshgrid(0:0.001:1, 0:0.001:1); %Para poder graficar bien la superficie
    
    z = f(xi, yi); %Función generadora del (0,0)
    
    triangular_region = (xi + yi <= 1); % Graficar solo el dominio de nuestra función
    
    z(~triangular_region) = NaN;
    
    figure;
    
    surf(xi, yi, z); %Grafico
    
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    colormap(jet);
    shading interp;
    alpha(0.7);  % Hacer la superficie semi-transparente
    
    % Añadir el plano z = 0 como otra superficie
    hold on; % Mantener el gráfico actual para agregar el plano
    
    z_plane = zeros(size(xi)); % Todos los puntos en z = 0
    
    plane = surf(xi, yi, z_plane, 'FaceColor', 'k', 'EdgeColor', 'none');
    alpha(plane, 0.3); % Hacer el plano semi-transparente para no bloquear la superficie
    
    grid on;      % Añadir la rejilla
    
    hold off
end