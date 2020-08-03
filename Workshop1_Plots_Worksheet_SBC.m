%% Matlab Series Workshop
% 6 Aug 2020 at SBC
% Data Visualization (Plots)
% Kan Kanjanapas (Ph.D)


clc;
close all;
clear all;


% 1. Plot

fs = 100;       % Sampling Frequency [Hz]
Ts = 1/fs;      % Sampling Time [s]

t_vec = [0:Ts:5]';  % Vector of time stamps

% sinusoidal waveform: A*sin(omega*t + phase_shift)

% Assume the reading from sensor A is given by:
f_1 = 1;    

x_1 = 1*sin(2*pi*f_1*t_vec + 0);
x_2 = 2*sin(2*pi*f_1*t_vec + 0);
x_3 = 3*sin(2*pi*f_1*t_vec + 0);
x_4 = 4*sin(2*pi*f_1*t_vec + 0);
x_5 = 5*sin(2*pi*f_1*t_vec + 0);

x_ceil = [];
for ii = 1:5
    x_ceil{ii} = ii*sin(2*pi*f_1*t_vec + 0);
end



%% 1) Plot

% Your first plot
close all

% Insert your code here **************

%%
% Refine the first plot 1.0

% Insert your code here **************


%%

% Refine the first plot: version 1.1 --------------------------------------------------------------------------------------------------
close all;
figure;
set( gcf, 'Position', [0 0 2560 1280]/2 );

% Insert your code here **************


xlabel('Time [s]');
ylabel('Position [m]');
title('Position Measurement from Sensor A');
grid on;
set(gca, 'FontSize', 14);




%%

% Refine the first plot: version 1.2 --------------------------------------------------------------------------------------------------
close all;
figure;
set( gcf, 'Position', [0 0 2560 1280]/2 );

% Line 1
% Insert your code here **************

hold on;

% Line 2
% Insert your code here **************

xlabel('Time [s]');
ylabel('Position [m]');
title('Position Measurement from Sensor A');

% Legend
% Insert your code here **************

grid on;
set(gca, 'FontSize', 14);



%%

% Refine the first plot: version 1.3 --------------------------------------------------------------------------------------------------
close all;
Color_Matrix = [0.0  0.0  1.0; 
                0.0  0.5  0.0;
                1.0  0.0  0.0;
                0.0  1.0  1.0;
                1.0  0.0  1.0];

MarkerEdgeColor_Ceil = {'r', 'k', 'b', 'm', 'g'};

figure;
set(gcf, 'Position', [0 0 2560 1280]/2);


for ii = 1:5
    
    % Insert your code here **************
    
    % Insert your code here ************** (if)
    
    
 
end
xlabel('Time [s]');
ylabel('Position [m]');
title(sprintf('Position Measurement from Sensor A with sampling frequency of %.1f [Hz]', f_1));
h_legend = legend('x_1', 'x_2', 'x_3', 'x_4', 'x_5'); % ************************************************************
set(h_legend, 'Location', 'NorthEast', 'Color', [1.0  1.0  0.9]);
grid on;
set(gca, 'FontSize', 14);



%%

% Refine the first plot: version 1.4 --------------------------------------------------------------------------------------------------

Color_Matrix = [0.0  0.0  1.0; 
                0.0  0.5  0.0;
                1.0  0.0  0.0;
                0.0  1.0  1.0;
                1.0  0.0  1.0];
            
MarkerEdgeColor_Ceil = {'r', 'k', 'b', 'm', 'g'};

close all;
figure;
set(gcf, 'Position', [0 0 2560 1280]/2);
for ii = 1:5
    plot(t_vec, x_ceil{ii}, 'LineStyle', '-', 'LineWidth', 2, 'Color', Color_Matrix(ii,:), ...
     'Marker', 'o', 'MarkerEdgeColor', MarkerEdgeColor_Ceil{ii}, 'MarkerFaceColor', 0.5*[1 1 1]);
    if (ii == 1)
        hold on;
    elseif (ii == 5)
        line([0 max(t_vec)], [0 0], 'LineStyle', '-', 'LineWidth', 2, 'Color', 'k');  % ****************
    end
end
xlabel('Time [s]');
ylabel('Position [m]');
title(sprintf('Position Measurement from Sensor A with sampling frequency of %.1f [Hz]', f_1));
h_legend = legend('x_1', 'x_2', 'x_3', 'x_4', 'x_5');
set(h_legend, 'Location', 'NorthEast', 'Color', [1.0  1.0  0.9]);
grid on;
set(gca, 'FontSize', 14);
axis([0 1 min(x_ceil{5}) max(x_ceil{5})]);


%% 2) Subplot

% Subplot Version 2.1 ----------------------------------------------------------------------------------------------------------
close all;
figure;
set(gcf, 'Position', [0 0 2560 1280]/2);

% Insert your code here **************


%%

% Subplot Version 2.2 ----------------------------------------------------------------------------------------------------------
% figure;
% set(gcf, 'Position', [0 0 2560 1280]/2);
% 
% subplot(2,2,1);
% plot(t_vec, x_1, 'r');
% 
% 
% subplot(2,2,2);
% plot(t_vec, x_2, 'g');
% 
% 
% subplot(2,2,3);
% plot(t_vec, x_3, 'b');
% 
% 
% subplot(2,2,4);
% plot(t_vec, x_4, 'k');


% Subplot Version 2.3 ----------------------------------------------------------------------------------------------------------
figure;
set(gcf, 'Position', [0 0 2560 1280]/2);

subplot(2,2,[1 2]);
plot(t_vec, x_1, 'r');


subplot(2,2,3);
plot(t_vec, x_2, 'g');



subplot(2,2,4);
plot(t_vec, x_4, 'k');





%% 3) Plot3


theta  = [0: pi/100: 20*pi]';
t_vec2 = [0 : 1 :length(theta)-1]'*Ts;
lambda = 0.2;

r = exp(-lambda*t_vec2);
x = r.*cos(theta);
y = r.*sin(theta);
z = theta/4;

% 3.1: Plot3 ----------------------------------------------------------------------------------------------------------
figure;
set(gcf, 'Position', [0 0 2560 1280]/2);
plot3(x, y, z, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
set(gca, 'FontSize', 16);
grid on;


%%
% 3.2: Plot3 Multiple View ----------------------------------------------------------------------------------------------------------

figure;
set(gcf, 'Position', [0 0 2560 1280]/2);

subplot(2,2,1);
plot3(x, y, z, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
set(gca, 'FontSize', 16);
grid on;


subplot(2,2,2);
plot3(x, y, z, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('X-Z Plane');
set(gca, 'FontSize', 16);
grid on;
view(0,0);


subplot(2,2,3);
plot3(x, y, z, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('Y-Z Plane');
set(gca, 'FontSize', 16);
grid on;
view(90,0);


subplot(2,2,4);
plot3(x, y, z, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('X-Y Plane');
set(gca, 'FontSize', 16);
grid on;
view(0,90);




%% 4) Surface

[X,Y] = meshgrid((-5:0.1:5)*pi, (-1:0.1:1)*pi);  % check size(X), size(Y)
Z = X.^2 + Y.^4;  % check size(Z)

% 4.1) Surf ----------------------------------------------------------------------------------------------------------
figure;
set(gcf, 'Position', [0 0 2560 1280]/2);

surf(X,Y,Z);
colorbar;

xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('Surf XYZ');
set(gca, 'FontSize', 14);
grid on;

% 4.2) Surf (Cont) ----------------------------------------------------------------------------------------------------------

% figure;
% set(gcf, 'Position', [0 0 2560 1280]/2);
% 
% surf(X,Y,Z, 'FaceAlpha', 0.5, 'EdgeColor', 'r');
% colorbar;
% 
% xlabel('X [-]');
% ylabel('Y [-]');
% zlabel('Z [-]');
% title('Surf XYZ');
% set(gca, 'FontSize', 14);
% grid on;



%% 5) Contour (contour3, countourc, contourf, quiver)

% [X,Y] = meshgrid((0:0.1:5)*pi, (0:0.1:1)*pi);  % check size(X), size(Y)
% Z = X.^2 + Y.^4;  % check size(Z)


figure;
set(gcf, 'Position', [0 0 2560 1280]/2);

subplot(2,2,1);
contour(X, Y, Z, 'ShowText', 'on');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('Contour');
set(gca, 'FontSize', 14);
grid on;

subplot(2,2,2);
contour3(X, Y, Z, 'ShowText', 'on');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('Contour3');
set(gca, 'FontSize', 14);
grid on;

subplot(2,2,3);
contourf(Z, 'ShowText', 'on');
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('Contourf');
set(gca, 'FontSize', 14);
grid on;


[DX, DY] = gradient(Z, 0.2, 0.2);
subplot(2,2,4);
contour(X,Y,Z);
hold on;
quiver(X,Y,DX,DY);
hold off;
xlabel('X [-]');
ylabel('Y [-]');
zlabel('Z [-]');
title('Quiver');
set(gca, 'FontSize', 14);
grid on;




%% 6) Semilogx, Semilogy, loglog

X = logspace(1,3,100);
Y = X.^2;

figure;
set(gcf, 'Position', [0 0 2560 1280]/2);

subplot(1,3,1);
semilogx(X,Y, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
title('Semilogx');
set(gca, 'FontSize', 14);
grid on;

subplot(1,3,2);
semilogy(X,Y, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
title('Semilogy');
set(gca, 'FontSize', 14);
grid on;

subplot(1,3,3);
loglog(X,Y, 'LineWidth', 2, 'Color', 'b');
xlabel('X [-]');
ylabel('Y [-]');
title('LogLog');
set(gca, 'FontSize', 14);
grid on;


