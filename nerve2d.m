%Morris-Lecar equations, fitted to the human node of Ranvier. 

istim = 0.015; %(mA/cm2)

gL = 0.001; %S/cm2
gNa = 0.06; %S/cm2

VR = -84; % Reversal potential for leak current and resting potential. Schwarz JR1, Reid G, Bostock H. Pflugers Arch. 1995 Jun;430(2):283-92.
cm = 2*10^-6;% Membrane capacitanc 2 uF/cm2. FRANKENHAEUSER B, HUXLEY AF. J Physiol. 1964 Jun;171:302-15.
eNa = 65; %Reversal potential for Na current, mV

tspan = 1;              % Time span in sec
dt = 0.00001;             % Time step for forward Euler method
loop  = ceil(tspan/dt);   % No. of iterations of Euler

V = zeros(loop,1);
h = zeros(loop,1);
V(1) = VR;
h(1) = hss(VR);
finish=ceil(.1/dt);

for i=1:loop-1
    V(i+1) = V(i) + dt*(gNa*mss(V(i))*h(i)*(eNa-(V(i))) + gL*(VR-(V(i)))+istim)/(cm);
    h(i+1) = h(i) + dt*((hss(V(i))-h(i))/(tauh(V(i))));
    if i>finish
        istim=0.0143;
    end

end

plot(1000*(0:dt:(tspan-dt)),V)
ylim([-85 25])

function m = mss(V)
mMid = -39; %mV
mSlope = 5; %mV
m = 0.5*(1+tanh((V-mMid)/(2*mSlope)));
end

function h = hss(V)
hMid = -79.1; %mV Schwarz JR1, Reid G, Bostock H. Pflugers Arch. 1995 Jun;430(2):283-92.
hSlope = 7.6; %mV Schwarz JR1, Reid G, Bostock H. Pflugers Arch. 1995 Jun;430(2):283-92.
h = 0.5*(1+tanh(-(V-hMid)/(2*hSlope)));
end

function tau = tauh(V)
hMid=-79.1; %mV
hSlope=7.6; %mV
kh=100;
tau = 1/(kh*cosh((V-hMid)/(2*hSlope)));
end
