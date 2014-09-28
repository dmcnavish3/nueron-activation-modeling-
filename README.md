nueron-activation-modeling-
===========================
totaltime = 100; %total time
numsteps = 10000; %total divisions of the total time
step1 = totaltime/numsteps; % step size
t = 0:step1:totaltime % vector of all time spaces %

%constants
    cond_K = 36; %mS/cm^2 
    cond_Na = 120; %mS/cm^2
    cond_leak = 0.3; %mS/cm^2
    Volt_K = -12; %mV
    Volt_Na = 115; %mV
    Volt_leak = 10.6; %mV
    Volt_resting = -70; %mV
    Cm = 1; %microFarads/cm^2
    Vm = Volt_resting;
    
    %initial conditions to variables that are changing 
    alpham = 0.1*(( 25-Vm)/(exp((25-Vm)/10)-1));
    betam = 4*exp(-Vm/18);
    alphan =0.01*( 10-Vm)/(exp((10- Vm)/10)-1);
    betan =.125 * exp(-Vm/80);
    alphah =0.07 * exp(-Vm/20);
    betah = 1/(exp(( 30-Vm)/10)+1);
    
    %initial m,h,n gate variables
    m = alpham/(alpham + betam);
    n = alphan/(alphan + betan);
    h = alphah/(alphah + betah);
    
    %pulse that I will give it to start the action potential
    I = 0;
    
    %this is needed to graph the conductance later on. Just make a vetor of
    %running conductances 
    gNa= m^3*h*cond_Na;
    gK= n^4*cond_K;
    
for t = 1:numsteps %for every step time, re-evaluate the model (0.01 ms) 
        
    %gating variables which change with changing Vm
    alpham = 0.1.*((25 - Vm(t))/(exp((25 -  Vm(t))/10)-1));
    betam = 4.*exp(-Vm(t)/18);
    alphan = 0.01.*(( Vm(t)-10)/(exp((10 - Vm(t))/10)-1));
    betan = .125.*exp(-Vm(t)/80);
    alphah = 0.07 * exp(-Vm(t)/20);
    betah =1/(exp((30 - Vm(t))/10)+1);
    
    %currents based on conductance and voltage
    INa = (m(t)^3).*cond_Na .* h(t).* (Vm(t) - Volt_Na);
    IK = (n(t)^4).*cond_K .* (Vm(t) - Volt_K);
    Ileak = cond_leak .* (Vm(t) - Volt_leak);
    
    Iion = I - IK - INa - Ileak;
    
    Vm(t+1) = Vm(t)+step1.*(Iion./Cm);    

    m(t+1) = (m(t)+step1*(alpham.*(1-m(t))-(betam*m(t))));
    n(t+1) = (n(t)+step1*(alphan.*(1-n(t))-(betan*n(t))));
    h(t+1) = (h(t)+step1*(alphah.*(1-h(t))-(betah*h(t))));
 
    
    %make a vector of these conductances to graph later on
    gNa = [gNa (m(t).^3).*cond_Na*h(t)];
    gK = [gK (n(t).^4).*cond_K];
end

figure 
        plot(0:step1:100, Vm)
        title('Membrane potential vs. time for pulse into a neuron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)') 
        
        figure
        plot(0:step1:100, (gNa./gK))
        title('gNa/gK vs. time pulse into a nueron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('coduction of Na/ conduction of Potassium')

    
    
