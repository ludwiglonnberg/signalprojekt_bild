% Ladda data-filer
load('rawData.mat');
load('hrf_1ms.mat'); 

% Centrera k-space
k_space = fftshift(k_space);

% Skapa spatiala bilder från 3D raw data och invers fourier-transform
iFFT_images = fftshift(ifftn(ifftshift(k_space)));

% Gaussiskt lågpassfilter för att filtrera brus (justera a?)
[x, y, z] = meshgrid(linspace(-1, 1, 128), linspace(-1, 1, 128), linspace(-1, 1, 37));
a = 0.1;
filter = exp(-(x.^2*a^2 + y.^2*a^2 + z.^2*a^2) / (2*a^3));

% Applicera filtret i frekvensplanet
k_space_filter = k_space .* filter;
iFFT_filtered = fftshift(ifftn(ifftshift(k_space_filter)));

% Aktiveringssignal
fs = 0.5; % Samplingsfrekvens
T = 320; % Total tid
N = 160; % Antal prover
time_vector = (0:N-1) * 2; % Tid i sekunder
activation_times = [0 32.4 64.8 97.2 129.7 162.1 194.5 226.9 259.3 291.7]; % Aktiveringstider
activation_signal = zeros(1, N);

% Skapar en aktiveringssignal som varar i 16s där värdet sätt till 1
for t = activation_times
    idx = find(time_vector >= t & time_vector < t + 16);
    activation_signal(idx) = 1;
end

% Skalar om och konvolvera aktiveringssignalen med HRF (beskriver den förväntade hjärnaktiviteten i tid)
hrf = hrfScaled(1,1:N);
activation_convolved = conv(activation_signal, hrf, 'same');

% Omforma fMRI data till en 2D matris så att voxlar blir en funktion av tiden
fMRI_matrix = reshape(abs(iFFT_filtered), [], 160);

% Hitta korrelationen mellan voxlar och aktiveringssignal
correlation_map = corr(fMRI_matrix', activation_convolved');

% Omforma korrelationsmappen till samma storlek som originalbilden
correlation_map = reshape(correlation_map, [128, 128, 37]);

% Ange tröskelvärde för korrelationen (bestäm gräns för aktiverade voxlar, justera threshold?)
threshold = 0.3;
activated_voxels = correlation_map > threshold;

% Skapa en meshgrid med (bild)storleken för varje dimension
[x_coords, y_coords, z_coords] = meshgrid(1:128, 1:128, 1:37);

% Skala om koordinaterna m.a.p voxelstorleken
x_coords = x_coords * 1.875;  
y_coords = y_coords * 1.875;
z_coords = z_coords * 3.9; 

% Plotta 3D-figur för aktiverade voxlar
figure;
isosurface(x_coords, y_coords, z_coords, activated_voxels, 0.5);  % 0.5 tröskelvärde (striktare, eller ska vi ha samma på båda?)
colormap jet;
colorbar;
title('Aktiverade voxlar markerade');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
axis equal;
view(3);



% Öka den spatiala upplösningen genom “zero-padding”. Använd Matlabfunktionen padarray. 
    % Spec_Im_padded=padarray(Spec_Im,[x y z]);

% Undersök correlationen bättre så vi kan redovisa det

% Visualisera mer resultat