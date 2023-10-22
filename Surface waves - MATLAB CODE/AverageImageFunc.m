function Iavg = AverageImageFunc(I)
    numberOfImages = length(I); %Number of images processed in batch processor
for k = 1 : numberOfImages

    thisImage = double(I{k});

    if k == 1
        sumImage = thisImage;
    else
        sumImage = sumImage + thisImage;
    end
end

Iavg = sumImage / numberOfImages;
Iavg = uint8(Iavg); %converts class from double back to uint8

end