%% Give f a value here so you dont have to do the calulations every time when testing. Delete this section when finished!
% Calculate the scaling factor, f, by uncommenting the section starting at
% line ~50 in this file
f = 3.0400e-04;
%% Declare variables 
%close all
%clc
%Choose an image to look at when displaying only one image
%image nr a is used to find origin
a = 100;
a = a-20;


H = 6e-2; %water depth [m]

%framerate
%fps = 506;
fps = 300;

dt = 1/fps;

%Decide how fast the animations will update 
n = 1/500;

%Exclude radii smaller than r_min
r_min = 100;

%start and end time for marbles
%for 1cm depth. use thres=23. and AverageImageFunc(I(:))
% t1=55; 
% t2=230; 

% %for 6cm.use thres=25. and AverageImageFunc(I(:))
% t1=20;
% t2=250;

% %for 10cm. thres=25. and AverageImageFunc(I(1:15))
% t1=20; 
% t2=170;

% %for 6cm offcenter. a=70. thres=25 and AverageImageFunc(I(:))
% t1=20;
% t2=250;

% %start and end time for droplets
% %for 1cm depth. a=70, thres=25 and AverageImageFunc(I(:))
%t1=30;
%t2=200;

% %for 3cm. a=70. thres=25 and AverageImageFunc(I(:))
% t1=1;
% t2=170;

% %for 6cm. a=70. thres=18 and Iavg = AverageImageFunc(I(185:200));%works for marble depth 1cm, 6cm
t1=30;
t2=130; 




g = 9.82; %Gravitational acceleration [m/s^2]
rho = 997; %Water density [kg/m^3]

%Surface tension, sigma, given by the equation from Vargatif "International
%tables of the surface tension of water". https://srd.nist.gov/jpcrdreprint/1.555688.pdf
%Constants given in the text by Vargaftik:
Tc = 647.15; %[K] 
B = 235.8e-3; %[N/m]
b = -0.625; %[-]
mu = 1.256; %[-]
%My measured water temperature in Kelvin, T:
T = 22+273.15;
sigma = B* ( ( (Tc-T)/Tc )^mu ) * ( 1 + b*(Tc-T)/Tc); % Surface tension [Pa]
%--------------------------------------------------------------------------

%% Calculate the scaling factor to convert measurements from pixels to meter
% % % 
% % %Load three reference-length-images
% ImgLref1 = imread('/Users/Albin/Desktop/examensarbete/Final experiment/Images/Experiment runs/calibration_final_exp_20120823_113502_AM/calibration_final_exp0000000001.tif');
% ImgLref2 = imread('/Users/Albin/Desktop/examensarbete/Final experiment/Images/Experiment runs/calibration_final_exp_20120823_113502_AM/calibration_final_exp0000000002.tif');
% ImgLref3 = imread('/Users/Albin/Desktop/examensarbete/Final experiment/Images/Experiment runs/calibration_final_exp_20120823_113502_AM/calibration_final_exp0000000003.tif');
% 
% %Assemble ref.images in cell
% ImgLref = {ImgLref1,ImgLref2,ImgLref3};
% 
% figure;
% disp('Measure the ruler in the image from 1-29cm three times to get reference length in pixels')
% for i = 1:3
%     imshow(ImgLref{i})
%     d = drawline;
%     pos = d.Position;
%     diffPos = diff(pos);
%     diameter(i,:) = hypot(diffPos(1),diffPos(2));
% end
% L_meters = 28e-2;%[m]
% L_pixels = mean(diameter); % reference length in [pixels]
% f = L_meters / L_pixels; %scaling factor [m/px]



