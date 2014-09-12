nueron-activation-modeling-
===========================
%constants
cond_K = 36; %mS/cm^2 
cond_Na = 120; %mS/cm^2
Cond_leak = 0.3; %mS/cm^2
Volt_K = -12; %mV
Volt_Na = 115; %mV
Volt_leak = 10.6; %mV
Volt_resting = -70; %mV

figure %constant Vm
    plot(0:0.5:100, -70*(ones(length(0:0.5:100))))
    title('Membrane potential vs. time for steady state neuron')
    axis([0,100,-100,100])
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    
figure %constant conductance
    plot(0:0.5:100, (cond_Na./ cond_K)*(ones(length(0:0.5:100))))
    title('gNa/gK vs. time for steady state nueron')
    axis([0,100,0,10])
    xlabel('Time (ms)')
    ylabel('coduction of Na/ confuction of Potassium')
    
for s = 1:length(1:0.5:100) % s is step time for the whole time
    
    
    
    figure %membrane potential after a step pulse 
        plot(0:0.5:100, Vm)
        title('Membrane potential vs. time for pulse into a neuron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
    
        figure
        plot(0:0.5:100, (cond_Na./cond_K
        title('gNa/gK vs. time pulse into a nueron')
        axis([0,100,0,10])
        xlabel('Time (ms)')
        ylabel('coduction of Na/ confuction of Potassium')
    
    
end
