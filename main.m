clear; close all;

%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm segments
P = 400; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis
axle_spacing = 100; % Distance between train axles [mm]
num_axles = 6; % Total number of train axles
P_axle = P / num_axles; % Load per axle

%% 1. SFD, BMD under train loading
train_positions = linspace(0, L + axle_spacing * (num_axles - 1), n + 1); % Train movement positions
SFDi = zeros(length(train_positions), n+1); % SFD for each train position
BMDi = zeros(length(train_positions), n+1); % BMD for each train position

% Calculate SFD and BMD for each train position
for i = 1:length(train_positions)
    % Current position of the train's first axle
    train_start = train_positions(i);
    
    % Determine which axles are on the bridge
    axles_on_bridge = ((train_start + (0:num_axles-1) * axle_spacing) >= 0) & ...
                      ((train_start + (0:num_axles-1) * axle_spacing) <= L);
    
    % Axle positions on the bridge
    axle_positions = train_start + (0:num_axles-1) * axle_spacing;
    axle_positions = axle_positions(axles_on_bridge);
    
    % Define load distribution along the bridge
    w = zeros(1, n+1);
    for j = 1:length(axle_positions)
        axle_index = find(x >= axle_positions(j), 1);
        if ~isempty(axle_index)
            w(axle_index) = w(axle_index) - P_axle; % Apply axle load
        end
    end
    
    % Calculate SFD by numerical integration
    SFDi(i, :) = cumtrapz(x, w);
    
    % Calculate BMD by numerical integration of SFD
    BMDi(i, :) = cumtrapz(x, SFDi(i, :));
end

% Calculate envelope for SFD and BMD
SFD = max(abs(SFDi), [], 1); % SFD envelope
BMD = max(abs(BMDi), [], 1); % BMD envelope

%% 2. Plot Results
figure;
subplot(2, 1, 1);
plot(x, SFD, 'b');
title('Shear Force Diagram (Envelope)');
xlabel('Position along the bridge (mm)');
ylabel('Shear Force (N)');
grid on;

subplot(2, 1, 2);
plot(x, BMD, 'r');
title('Bending Moment Diagram (Envelope)');
xlabel('Position along the bridge (mm)');
ylabel('Bending Moment (N·mm)');
grid on;

%% 2. Define Bridge Parameters
param = [0, 100, 1.27; 400, 100, 1.27; 800, 100, 1.27; L, 100, 1.27];
bft = interp1(param(:,1), param(:,2), x); % Top Flange Width
tft = interp1(param(:,1), param(:,3), x); % Top Flange Thickness

%% 3. Calculate Sectional Properties
ybar = 50; % Centroidal axis height
I = bft .* tft.^3 / 12; % Moment of inertia (simplified)
Qcent = bft .* ybar; % First moment of area at centroid
Qglue = Qcent / 2; % First moment at glue joint

%% 4. Calculate Applied Stress
S_top = BMD .* ybar ./ I;
S_bot = -S_top;
T_cent = SFD .* Qcent ./ I;
T_glue = SFD .* Qglue ./ I;

%% 5. Material and Thin Plate Buckling Capacities
E = 4000; % Elastic modulus [MPa]
mu = 0.2; % Poisson's ratio
S_tens = 10; % Tensile strength [MPa]
S_comp = 8; % Compressive strength [MPa]
T_max = 5; % Shear strength [MPa]
T_gmax = 3; % Glue shear strength [MPa]

% Buckling capacities (arbitrary values for illustration)
S_buck1 = 0.8 * E * tft.^2 ./ (12 * (1 - mu^2));
S_buck2 = 0.6 * E * tft.^2 ./ (12 * (1 - mu^2));
S_buck3 = 0.4 * E * tft.^2 ./ (12 * (1 - mu^2));
T_buck = 0.5 * E * tft ./ (1 + mu);

%% 6. Factor of Safety (FOS)
FOS_tens = S_tens ./ abs(S_top);
FOS_comp = S_comp ./ abs(S_bot);
FOS_shear = T_max ./ abs(T_cent);
FOS_glue = T_gmax ./ abs(T_glue);
FOS_buck1 = S_buck1 ./ abs(S_top);
FOS_buck2 = S_buck2 ./ abs(S_top);
FOS_buck3 = S_buck3 ./ abs(S_top);
FOS_buckV = T_buck ./ abs(T_cent);

%% 7. Min FOS and the failure load Pfail
minFOS = min([FOS_tens; FOS_comp; FOS_shear; FOS_glue; FOS_buck1; FOS_buck2; FOS_buck3; FOS_buckV], [], 1);
Pf = P * min(minFOS);

%% 8. Vfail and Mfail
Mf_tens = S_tens .* I ./ ybar;
Mf_comp = S_comp .* I ./ ybar;
Vf_shear = T_max .* I ./ Qcent;
Vf_glue = T_gmax .* I ./ Qglue;
Mf_buck1 = S_buck1 .* I ./ ybar;
Mf_buck2 = S_buck2 .* I ./ ybar;
Mf_buck3 = S_buck3 .* I ./ ybar;
Vf_buckV = T_buck .* I ./ Qcent;

%% 9. Output plots of Vfail and Mfail
figure;
subplot(2,3,1);
hold on; grid on;
plot(x, Vf_shear, 'r');
plot(x, -Vf_shear .* SFD, 'k');
legend('Matboard Shear Failure');
xlabel('Distance along bridge (mm)');
ylabel('Shear Force (N)');

subplot(2,3,2);
hold on; grid on;
plot(x, Vf_glue, 'r');
plot(x, -Vf_glue .* SFD, 'k');
legend('Glue Shear Failure');
xlabel('Distance along bridge (mm)');
ylabel('Shear Force (N)');

subplot(2,3,3);
hold on; grid on;
plot(x, Vf_buckV, 'r');
plot(x, -Vf_buckV .* SFD, 'k');
legend('Matboard Shear Buckling Failure');
xlabel('Distance along bridge (mm)');
ylabel('Shear Force (N)');

subplot(2,3,4);
hold on; grid on;
plot(x, Mf_tens, 'r');
plot(x, -Mf_tens .* BMD, 'k');
legend('Matboard Tensile Failure');
xlabel('Distance along bridge (mm)');
ylabel('Bending Moment (Nmm)');

subplot(2,3,5);
hold on; grid on;
plot(x, Mf_buck1, 'r');
plot(x, -Mf_buck1 .* BMD, 'k');
legend('Matboard Buckling Failure, Top Flange');
xlabel('Distance along bridge (mm)');
ylabel('Bending Moment (Nmm)');

subplot(2,3,6);
hold on; grid on;
plot(x, Mf_buck3, 'r');
plot(x, -Mf_buck3 .* BMD, 'k');
legend('Matboard Buckling Failure, Web');
xlabel('Distance along bridge (mm)');
ylabel('Bending Moment (Nmm)');
