

clear, clc
% for NR = [1 3 6 10] %For marble
%for NR = [1 3 6] %For droplet
disp('OBS! The Matlab package "Image processing and computer vision" must be installed for the script to work')
disp('Remember to add folders to path "Images from experiment", and "data" in order to run this')
disp('Remember to change the filename ending to correct water depth')
disp('This is done in the matfile Calculations_alt0.m')
disp(' ')
disp('Also, remember to change to the correct depth H in declaringvariables.m')
disp('Also, remember to change the counter NR in the for loop below')
disp('STEPLENGTH FOR MEAN VALUES ARE 0.1 CM!!!!!!!!')
for NR = [6]


%% Loading variables
DeclaringVariables;
%a=120-50; %for marble in 10cm and 10cmADJUSTED water depth.
a=120; %for marble in 1 and 6cm water depth
%% Loading images
%I = load( sprintf('Images_MarbleOFFCENTER_%dcm.mat',NR) ).I; %Load the saved batch-processed images
I = load( sprintf('Images_Marble_%dcm.mat',NR) ).I; 
%I = load( sprintf('Images_MarbleADJUSTED_%dcm.mat',NR) ).I; 
%I = load( sprintf('Images_Droplet_%dcm.mat',NR) ).I; 
%%

% Create an averaged image-------------------------------------------------
%Iavg = AverageImageFunc(I(1:15)); %works for marble 10cm,
Iavg = AverageImageFunc(I(:));%works for marble depth 1cm, 6cm
%Iavg = AverageImageFunc(I(185:200));%works for marble depth 1cm, 6cm
%Iavg=AverageImageFunc(I(45:end));
%Threshold level for binarizing--------------------------------------------
%thres = 18; %for droplet with depth 6cm
%thres = 23; %for marble with depth 1cm
thres = 25; %for marble with depth 6cm


%Remove background, adjust contrast, threshold, edge-detection-------------
I2 = cell(1,height(I));
for i = 1:height(I)
    %Remove background from images-----------------------------------------
    I2{i} = imsubtract(I{i},Iavg) ;

    %Adjust image contrast-------------------------------------------------
    I3{i} = imadjust(I2{i});

    %Filter image----------------------------------------------------------
    I3{i} = imdiffusefilt(I3{i});
   
    %Threshold image-------------------------------------------------------
    I3{i} = I3{i} > thres;

    %Clean up image--------------------------------------------------------
    I3{i} = bwareaopen(I3{i},50) ;
    I3{i} = imclearborder(I3{i});

end

%% Calculate center of impact (origin for wave pattern)
[center,radius] = imfindcircles(I3{a},[125 320],'ObjectPolarity','bright', ...
    'Sensitivity',0.85, 'Method','twostage'); 

x0 = round(center(1));
y0 = round(center(2));


%% Run both the horizontal and vertical case and save the concatenated data
Calculations_alt0;

%% Plot the results (lambda vs c, and errors) without outliers

h = NR;
%plotResults_12Aug;
updated_plots_results;

% 
% matFilename_celldata = sprintf('celldata.mat');%Horizontal
% matFilename_celldata2 = sprintf('celldata2.mat'); %Vertical
% %Horizontal
% lam_i = load(matFilename_celldata).lam_i;
% c_exp_lead = load(matFilename_celldata).c_exp_lead;
% c_exp_trail = load(matFilename_celldata).c_exp_trail;
% 
% lam_i2 = load(matFilename_celldata2).lam_i;
% c_exp_lead2 = load(matFilename_celldata2).c_exp_lead;
% c_exp_trail2 = load(matFilename_celldata2).c_exp_trail;
% 
% testdata = [lam_i', c_exp_lead', c_exp_trail'];
% testdata2 = [lam_i2', c_exp_lead2', c_exp_trail2'];
% 
% % roundtestdata = cellfun(@(x) round(x,2,'significant'),testdata,'UniformOutput',false);
% % roundtestdata2 = cellfun(@(x) round(x,2,'significant'),testdata2,'UniformOutput',false);
% %Try skipping the rounding of data here:
% roundtestdata = testdata;
% roundtestdata2 = testdata2;
% 
% 
% bw = I3(t1:t2); 
% 
% FIGUR = figure;
% set(FIGUR, 'Color','w', 'Position', [0 300 1800 400])
% 
% %Keep legend from updating
% set(FIGUR,'defaultLegendAutoUpdate','off');
% 
% tlo = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
% nexttile
% 
% %for i =1:(t2-t1)
%     %for i =1:(t2-t1-30)
%     for i=1
%     subplot(1,2,1)
%     imshow(bw{i})
%     title(['Image nr:' num2str(i)])
%     hold on
%    % plot(xcrestrhs{i},ones(length(xcrestrhs{i}),1)*y0,'xr','LineWidth',2)
%    % plot(xcrestlhs{i},ones(length(xcrestlhs{i}),1)*y0,'xr','LineWidth',2)
%   
%    %x-dir
%    plot(x_horizontal_rhs{i+t1-1},ones(length(x_horizontal_rhs{i+t1-1}),1)*y0,'xr','LineWidth',2)
%    plot(x_horizontal_lhs{i+t1-1},ones(length(x_horizontal_lhs{i+t1-1}),1)*y0,'xr','LineWidth',2)
% 
%     %y-dir
%     plot(ones(length(y_vertical_rhs{i+t1-1}),1)*x0,y_vertical_rhs{i+t1-1},'xr','LineWidth',2)
%     plot(ones(length(y_vertical_lhs{i+t1-1}),1)*x0,y_vertical_lhs{i+t1-1},'xr','LineWidth',2)
% 
% 
%     subplot(1,2,2)
%     hold on
%     %Horizontal
%    plot(roundtestdata{i,1}*1e2,roundtestdata{i,2}*1e2,'.','MarkerSize',10,'Color','#0072BD')
%    plot(roundtestdata{i,1}*1e2,roundtestdata{i,3}*1e2,'.','MarkerSize',10,'color',	'#0072BD') % different color: '#D95319'
%     
%     %Vertical
%    plot(roundtestdata2{i,1}*1e2,roundtestdata2{i,2}*1e2,'.','MarkerSize',10,'Color','#D95319')
%    plot(roundtestdata2{i,1}*1e2,roundtestdata2{i,3}*1e2,'.','MarkerSize',10,'color','#D95319') % different color: '#D95319'
%     
%     %theoretical values
%     plot(lambda*1e2,c*1e2,'Color','#EDB120')
%     plot(lambda*1e2,c_cap*1e2,'Color','#77AC30')
%     plot(lambda*1e2,c_grav*1e2,'Color','#A2142F')
% 
%     pause(n+0.1);
%     legend('$c_{exp,horizontal}$','$c_{exp,vertical}$','$c_{mixed}$','$c_{capillary}$','$c_{gravity}$','Interpreter','latex')
%     title('Wavelength vs wavespeed','fontsize',15,'Interpreter','latex')
%     ylabel('$c$ [cm/s]','Interpreter','latex')
%     xlabel('$\lambda$ [cm]','Interpreter','latex')
%     xlim([0 6])
%     ylim([0 100])
%     set(gca,'fontsize',15)
%     
% 
% end
% %  
% 
% 
% 
% 






%clear, clc
end


%% Images below, can be deleted.

figure;
for i = 90:140
    imshow(I3{i})
    hold on

    viscircles(center,r_min);
    plot(x0,y0,'gx','LineWidth',2);
    %x-dir
    plot(x_horizontal_rhs{i},ones(length(x_horizontal_rhs{i}),1)*y0,'xr','LineWidth',2)
   plot(x_horizontal_lhs{i},ones(length(x_horizontal_lhs{i}),1)*y0,'xr','LineWidth',2)

    %y-dir
    plot(ones(length(y_vertical_rhs{i}),1)*x0,y_vertical_rhs{i},'xr','LineWidth',2)
    plot(ones(length(y_vertical_lhs{i}),1)*x0,y_vertical_lhs{i},'xr','LineWidth',2)
    pause(n+0.1);
    title(['Binary filtered image I3 with highlighted wavecrests for image nr:' num2str(i)])
    xlabel('r_i+x_0')
    ylabel('y_0')
    drawnow
end

% 
% 
% %%
% 
% figure; imshowpair(I{a},I2{a},'montage')
% title('Image nr a (left). Background subtracted image nr a (right)')
% 
% figure; imshow(I{a})
% title('Original image')
% 
% figure; imshow(I2{a})
% title('Background subtracted image')
% 
% figure; imshow(Iavg)
% title('Averaged image')
% % 
% %Display image with drawn circle
% figure; imshow(I3{t2})
% title('Background subtracted image with imfindcircles used to highlight origin.')
% hold on
% viscircles(center,radius);
% %Display origin of wave pattern
% plot(x0,y0,'yx','LineWidth',2);
% hold off
% %%
% 
% tlo = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% imshow(Iavg)
% nexttile
% imshow(I{a})
% nexttile
% imshow(I2{a})
% %%
% %figure;
% %montage({Iavg,I{a}, imadjust(I2{a})},'Size',[1 3])
% 
% 
% [J,rect] = imcrop(I{a}); 
% J2a=imcrop(imadjust(I2{a}),rect);
% Javg=imcrop(Iavg,rect);
% 
% %%
% 
%  FIGUR33 = figure;
%     set(FIGUR33, 'Color','w', 'Position', [0 300 1200 400])
% %montage({Iavg,I{a}, imadjust(I2{a})},'Size',[1 3])
% montage({Iavg,I{a}, I2{a}},'Size',[1 3])
% %            figName = sprintf('Iavg_I_I2.png');
% %   exportgraphics(FIGUR33,figName,'Resolution',800)
% 
% 
% 
% 
% 
% % crop and save I vs I3, first image in the results part
%     xtickvalues_mm=[5:5:25]*10;%[mm]
%     ytickvalues_mm=[5:5:25]*10;%[mm]
%     xtickvalues=xtickvalues_mm/(f*1000);
%     ytickvalues=ytickvalues_mm/(f*1000);
% 
%     xticks(xtickvalues)
%     xticklabels(string(xtickvalues_mm))
%     yticks(ytickvalues)
% 
%     yticklabels(string(ytickvalues_mm))
% 
% 
% [J,rect] = imcrop(I{a});
% 
% % I_i=imcrop(I{a},rect);
% I_i=imcrop(imadjust(I{a}),rect); 
% I3_i=imcrop(I3{a},rect);
% 
% figure;
% tlo = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
% nexttile
% imshow(I_i)
% xlabel('x [mm]','Interpreter','latex')
% ylabel('y [mm]','Interpreter','latex')
% axis on
% xticks(xtickvalues)
% xticklabels(string(xtickvalues_mm))
% yticks(ytickvalues)
% yticklabels(string(ytickvalues_mm))
% set(gca,'fontsize',15)
% 
% %----
% 
% nexttile
% imshow(I3_i)
% xlabel('x [mm]','Interpreter','latex')
% ylabel('y [mm]','Interpreter','latex')
% axis on
% xticks(xtickvalues)
% xticklabels(string(xtickvalues_mm))
% yticks(ytickvalues)
% yticklabels(string(ytickvalues_mm))
% set(gca,'fontsize',15)
% % 
% % figName = sprintf('IandI3_droplet_6cm.png');
% % exportgraphics(tlo,figName,'Resolution',800)
