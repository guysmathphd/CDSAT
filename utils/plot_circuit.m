function plot_circuit(A)

n = size(A,1);
t = 0:2*pi/n:2*pi-2*pi/n;
t = t';
xy = [cos(t),sin(t)];
figure;
gplot(A,xy,'o-');
end