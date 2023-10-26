rp = 6371+400;
u = 3.986004418e5;
vcs = sqrt(u/rp);
v = vcs + 3;

dv = 0.002;

h = rp*v;
p = h^2/u;
e = (rp*v^2)/u -1;
a = p/(1-e^2);

t = 0;

while a > rp
    v = v - dv;
    h = rp*v;
    p = h^2/u;
    e = (rp*v^2)/u -1;
    a = p/(1-e^2);

    TP = (2*pi/sqrt(u)) * a^(3/2);
    t = t + TP;
end

t / (365*24*3600)
