%



for i = h
    %Declare water depth
    H = i*1e-2;

    %Load data from marbles
       matFilename_horizontal = sprintf('DATA22_horizontal_Marble_%dcm.mat',i);
        matFilename_vertical = sprintf('DATA22_vertical_Marble_%dcm.mat',i);

       % matFilename_horizontal = sprintf('DATA22_horizontal_MarbleOFFCENTER_%dcm.mat',i);
        %matFilename_vertical = sprintf('DATA22_vertical_MarbleOFFCENTER_%dcm.mat',i);
   % matFilename_horizontal = sprintf('DATA22_horizontal_Droplet_%dcm.mat',i);
   % matFilename_vertical = sprintf('DATA22_vertical_Droplet_%dcm.mat',i);

    EXPDATA_horizontal = load(matFilename_horizontal);
    EXPDATA_vertical = load(matFilename_vertical);

    %Reduce number of significant numbers
    %     ExpdataTEST = { round(EXPDATA_horizontal.EXPDATA,2,'significant'), round(EXPDATA_vertical.EXPDATA,2,'significant')};
    ExpdataTEST = {EXPDATA_horizontal.EXPDATA, EXPDATA_vertical.EXPDATA};
    ExpdataTEST = [ExpdataTEST{1};ExpdataTEST{2}];
    % str = {'Horizontal','Vertical'};
    FIGUR = figure;
    set(FIGUR, 'Color','w', 'Position', [0 300 1200 400])

    EXPDATA = ExpdataTEST;


    %Calculaing the theoretical values
    lambda = linspace(0,0.15,1000);
    c = sqrt( ( g* lambda /(2* pi) + 2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
    c_cap = sqrt( (  2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
    c_grav = sqrt( ( g* lambda /(2* pi) ) .*tanh( 2*pi* H./lambda ) );

    % Find and remove outliers-----------------------------------------
    LAMVALS = [EXPDATA(:,1);EXPDATA(:,1)]; %Store wavelengths on top of each other
    CVALS = [EXPDATA(:,2);EXPDATA(:,3)];
    newDATA = [LAMVALS,CVALS]; %Put together matrix of trailing and leading crest wavelength and speed
    OutlRemovedDATA = rmoutliers(newDATA,'grubbs'); %Find outliers

    %         %------------------------------------------------------------------

    %Data without outliers---------------------------------------------
    Lamexp = OutlRemovedDATA(:,1)*1e-2; %multiply by 1e-2 to convert to [m] instead of [cm]
    Cexp = OutlRemovedDATA(:,2)*1e-2; %multiply by 1e-2 to convert to [m/s] instead of [cm/s]

    %Plot experimental data without outliers-----------------------------------
    plot(Lamexp*1e2,Cexp*1e2,'.','MarkerSize',10)
    hold on
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------Try to find mean values to include in plot here!!------------------
%--------------------------------------------------------------------------
lam_cm = [Lamexp]*1e2;%[cm] all lambda values in cm
c_cm=Cexp*1e2; %[cm/s] all c values in cm/s

steps = 0:0.1:round(max(lam_cm))+0.5;
for k=1:length(steps)
    lam_cm1=steps(k);
    lam_cm2 = lam_cm1+0.1;
    c_average(k,:) = [(lam_cm1+lam_cm2)/2, mean(c_cm((lam_cm>=lam_cm1) & (lam_cm<lam_cm2)))];

end
plot(c_average(:,1),c_average(:,2),'.r','MarkerSize',17)
hold on
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
 title(['Wavelength vs wavespeed without outliers. Waterdepth = ',num2str(i),'cm'],'fontsize',17,'Interpreter','latex')
   % title(['Wavelength vs wavespeed without outliers. Waterdepth = 6cm'],'fontsize',15,'Interpreter','latex')
    
   ylabel('$c$ [cm/s]','Interpreter','latex')
    xlabel('$\lambda$ [cm]','Interpreter','latex')
    xlim([0 6])
    ylim([0 100])

    


    %Plot theoretical values of lambda vs dR/dt--------------------------------
    plot(lambda*1e2,c*1e2,lambda*1e2,c_cap*1e2,lambda*1e2,c_grav*1e2,'linewidth',1.5)
    legend('$c_{experimental}$','$c_{mean}$','$c_{mixed}$','$c_{capillary}$','$c_{gravity}$','Interpreter','latex')
   set(gca,'fontsize',17)
    hold off
   
    %Save the figures as png files
    figName = sprintf('Results_Droplet_%dcm.png',i);

         %  figName = sprintf('Results_Marble_%dcm.png',i);
%figName = sprintf('Results_MarbleOFFCENTER_%dcm.png',i);
 %exportgraphics(FIGUR,figName,'Resolution',800)




    FIGUR2 = figure;
    set(FIGUR2, 'Color','w', 'Position', [0 300 1200 400])

    %Calculate the theoretical wavespeeds
    % scatter absolute error vs lambda and include trendline  (|c_exp-c_theory|)
    %Theoretical wavespeeds----------------------------------------------------
    c_theory = @(x) sqrt( ( g* x /(2* pi) + 2*pi*sigma./ (rho *x) ) .*tanh( 2*pi* H./x ) ); %[m/s]
    c_cap_theory = @(x) sqrt( (  2*pi*sigma./ (rho *x) ) .*tanh( 2*pi* H./x ) ); %[m/s]
    c_grav_theory = @(x) sqrt( ( g* x/(2* pi) ) .*tanh( 2*pi* H./x ) ); %[m/s]

    % %Absolute normalized errors:--------------------------------------------------------
    err_cap = abs( Cexp - c_cap_theory( Lamexp ) )./c_cap_theory( Lamexp ) ;
    err_grav = abs( Cexp - c_grav_theory( Lamexp ) )./c_grav_theory( Lamexp ) ;
    err_capgrav = abs( Cexp - c_theory( Lamexp ) )./ c_theory( Lamexp );

    %scatter-plot lambda vs relative absolute errors
    %         subplot(1,2,j)
    plot(Lamexp*1e2,err_capgrav,'.','MarkerSize',8)
    hold on
    plot(Lamexp*1e2,err_cap,'.','MarkerSize',8)
    plot(Lamexp*1e2,err_grav,'.','MarkerSize',8)

    %plot trendlines-----------------------------------------------------------
    [px,py] = trendline(Lamexp,err_capgrav);
    plot(px*1e2,py,'Color','#0072BD','linewidth',1.5)

    [px,py] = trendline(Lamexp,err_cap);
    plot(px*1e2,py,'color',	'#D95319','linewidth',1.5)

    [px,py] = trendline(Lamexp,err_grav);
    plot(px*1e2,py,'color','#EDB120','linewidth',1.5)



    xlabel('$\lambda$ [cm]','Interpreter','latex')
   % ylabel('$|c_{exp}-c_{theory}|/c_{theory}$ [-]','Interpreter','latex')
    ylabel('$\epsilon$ [-]','Interpreter','latex')
    
    legend('Capillary-Gravity','Capillary','Gravity','','','','Interpreter','latex')
   title(['Absolute normalized errors vs wavelength with trendlines. Waterdepth = ',num2str(i),'cm'],'Interpreter','latex')
% title(['Absolute normalized errors vs wavelength with trendlines. Waterdepth = 6cm'],'Interpreter','latex')

    set(gca,'fontsize',17)

    ylim([0 5])
    xlim([0 6])
    hold off


%     %     %Save the figures as png files
figName2 = sprintf('Errors_Droplet_%dcm.png',i);
         %figName2 = sprintf('Errors_Marble_%dcm.png',i);
%figName2 = sprintf('Errors_MarbleOFFCENTER_%dcm.png',i);
%exportgraphics(FIGUR2,figName2,'Resolution',800)


end



