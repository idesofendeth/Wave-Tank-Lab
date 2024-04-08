function phaseSpeed = phaseSpeedCalc(Radius1,Radius2,dt)
%phaseSpeedCalc: Finds the phase speed based on two radius vectors from two
%sequential images. Radius2 is the radii from image time t+1 and Radius1
%is the radii from time t
%   Detailed explanation goes here

diffr=[];
phaseSpeed=[];
j=1; %index for phase speeds
for i=1:length(Radius2)
    diffr=Radius2(i)-Radius1
    if (any(diffr>0))
        [val,idx]=min(abs(diffr))
        phaseSpeed(j)=(Radius2(i)-Radius1(idx))/dt
        j=j+1;
    end

end
phaseSpeed=phaseSpeed';
end