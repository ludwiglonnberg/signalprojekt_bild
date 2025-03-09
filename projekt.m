% Ladda data-filer
load('rawData.mat');
load('hrf_1ms.mat'); 

% Centrera k-space
kSpace = fftshift(kSpace);

% Skapa spatiala bilder från 3D raw data och invers fourier-transform
iFFT_images = fftshift(ifftn(ifftshift(kSpace)));

% Gaussiskt lågpassfilter för att filtrera brus (justera a?)
[x, y, z] = meshgrid(linspace(-1, 1, 128), linspace(-1, 1, 128), linspace(-1, 1, 37));
a = 0.08;
filter = exp(-(x.^2*a^2 + y.^2*a^2 + z.^2*a^2) / (2*a^3));

% Applicera filtret i frekvensplanet
k_space_filter = kSpace .* filter;
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
correlation_values = corr(fMRI_matrix', activation_convolved');

% Omforma korrelationsmappen till samma storlek som originalbilden
correlation_map = reshape(correlation_values, [128, 128, 37]);

% Ange tröskelvärde för korrelationen (bestäm gräns för aktiverade voxlar, justera threshold?)
threshold = 0.5;
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

% Öka den spatiala upplösningen genom “zero-padding”. Använd
% Matlabfunktionen padarray. (Detta verkar vara VG-nivå, ska vi göra det?)
    % Spec_Im_padded=padarray(Spec_Im,[x y z]);


% Undersök correlationen bättre så vi kan redovisa det
% (känns som att detta är fett lågt, har jag gjort fel?
mean_corr = mean(correlation_values);
std_corr = std(correlation_values);
disp(['Medelvärde: ', num2str(mean_corr)]);
disp(['Standardavvikelse: ', num2str(std_corr)]);

% Visualisera mer resultat

%plottar mittenslicen av ofiltrerad och filtrerad
mid_slice = round(size(iFFT_images,3)/2);
figure;
subplot(1,2,1);
imagesc(abs(iFFT_images(:,:,mid_slice))); 
title('Ofiltrerad');

subplot(1,2,2);
imagesc(abs(iFFT_filtered(:,:,mid_slice))); 
title('Filtrerad');

figure;
histogram(correlation_values,50); 
xlabel('Korrelationsvärde');
ylabel('Antal voxlar');
title('Korrelation mellan voxlar och aktiveringssignal');

%visar alla filtrerade slices i z-axeln
for slice_idx = 1:37
    subplot(6, 7, slice_idx);
    imagesc(abs(iFFT_filtered(:,:,slice_idx))); 
    axis equal; axis off;
    title(['Slice ', num2str(slice_idx)]);
  
end    
