% testing how find peaks deals with two points at the same height

x=1:6;
y=[1 2 3 3 2 1];


[peaks locs]=findpeaks(y);

figure;
hold on
plot(x,y)
plot(locs,peaks,'*')