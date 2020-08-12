% An example of 3-D plot
% surf_eg.m
x = linspace(-1,1);
y = linspace(0,1);
[X,Y] = meshgrid(x,y);
Z = X.^2 + Y;
surf(X,Y,Z)

% plot3_eg.m
t = linspace(0,10*pi);
plot3(sin(t),cos(t),t);
grid on
axis square
% create a movie
fps = 15;

t = linspace(0,6,6*fps);
x = linspace(0,2*pi);
h = figure(1);
clear M;
% plot and save the frames
for i = 1:length(t)
    plot(x, sin(pi*t(i))*sin(x));
    axis([0,2*pi,-1,1]);
    M(i) = getframe(h);
end
movie2avi(M, ?movie_eg.avi?, ?compression?,?none?, ?fps?,fps);

