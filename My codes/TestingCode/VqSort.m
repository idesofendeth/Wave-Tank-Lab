function [XPeakLoc,YPeakLoc] = VqSort(Vqpeaks, Vqlocs, Vqwidths, Vqproms,VqNorm,theta, spaceThresh,R)
%VqSort Organizes and finds the locations of the peaks in the data.
% input: output from findpeaks
% output: fixed version of locations of peaks

%index initialization
k=2;
i=1;



while k<=length(Vqlocs)
    %diff of location positions
    Vqlocdiff=Vqlocs(k)-Vqlocs(k-1);
    % windows of location positions and peak positions one is comparing
    VqLocWindow=[Vqlocs(k-1) Vqlocs(k)];
    VqPeakWindow=[Vqpeaks(k-1) Vqpeaks(k)];
    VqwidthsWindow=[Vqwidths(k-1) Vqwidths(k)];
    VqpromsWindow=[Vqproms(k-1) Vqproms(k)];
    % filtering out locations of peaks that are too close to eachother
    % taking the largest peak to save, throwing away the other
    % there is some bug here, some of the points arent disappearing..
    if abs(Vqlocdiff)<spaceThresh
        %finding the value and index of the maximum between the
        %investigated points
        [VqPeakFilter I]=max(VqPeakWindow);
        %edge cases for start and end of vector
        if k==2
            VqlocsFixed=[VqLocWindow(I) Vqlocs(k+1:end)];
            VqPeaksFixed=[VqPeakFilter Vqpeaks(k+1:end)];
            VqwidthsFixed=[VqwidthsWindow(I) Vqwidths(k+1:end)];
            VqpromsFixed=[VqpromsWindow(I) Vqproms(k+1:end)];
        elseif k==length(Vqlocs)
            VqlocsFixed=[Vqlocs(1:k-2) VqLocWindow(I)];
            VqPeaksFixed=[Vqpeaks(1:k-2) VqPeakFilter];
            VqwidthsFixed=[Vqwidths(1:k-2) VqwidthsWindow(I)];
            VqpromsFixed=[Vqproms(1:k-2) VqpromsWindow(I)];
        else %otherwise do this
            VqlocsFixed=[Vqlocs(1:k-2) VqLocWindow(I) Vqlocs(k+1:end)];
            VqPeaksFixed=[Vqpeaks(1:k-2) VqPeakFilter Vqpeaks(k+1:end)];
            VqwidthsFixed=[Vqwidths(1:k-2) VqwidthsWindow(I) Vqwidths(k+1:end)];
            VqpromsFixed=[Vqproms(1:k-2) VqpromsWindow(I) Vqproms(k+1:end)];

        end
    else
        VqlocsFixed=Vqlocs;
        VqPeaksFixed=Vqpeaks;
        VqwidthsFixed=Vqwidths;
        VqpromsFixed=Vqproms;


    end
    Vqlocs=VqlocsFixed;
    Vqpeaks=VqPeaksFixed;
    Vqwidths=VqwidthsFixed;
    Vqproms=VqpromsFixed;

    k=k+1;

    %end of while loop
end
%fitting a polynomial to the data
% troubleshooting output!
% Vqlocs
% VqlocsFixed
% Vqproms
% VqpromsFixed
% Vqwidths
% VqwidthsFixed
%
PolyPeakLocs=[];
%troubleshooting figure
% figure(97)
% plot(R,VqNorm,'LineWidth',2)
% plot(VqlocsFixed,VqPeaksFixed,'*')

for ii=1:length(Vqlocs) %PREVIOUSLY VqlocsFixed, change if this fucks up, same with widths
    [fitresult, gof] = createPoly2FitV3(R, VqNorm,Vqlocs(ii),abs(Vqwidths(ii))/2);

    PolyPeakLocs=[PolyPeakLocs fitresult.b];

end
VqlocsFixed=PolyPeakLocs;

% Extract x and y positions from the peak radial locations
XPeakLoc=VqlocsFixed*cosd(theta);
YPeakLoc=VqlocsFixed*sind(theta);
% add values to x and y peak vectors


%troubleshooting polyfit plot
figure(98),clf
hold on
a=length(R(1:end))
plot(VqlocsFixed,VqPeaksFixed,'*')
plot(R,VqNorm,'LineWidth',2)
xlim([0.05 0.15])
%plot(R,VqNorm,'g--','LineWidth',2)
hold off
% %pause

%end of for loop
end
