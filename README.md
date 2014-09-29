nueron-activation-modeling-
===========================
%% part 1
clear
%set the time that this function will run for
totaltime = 100; 
%set the step of each re-evaluation of the derivative functions
step1 = .01;
%make a vector of the times that this function will be re-evaluated (in
%milliseconds) starting at zero. I will use this later to plot my points
%and to run the loops.
time = 0:step1:totaltime;

%constants that were given
cond_K = 36; %mS/cm^2 
cond_Na = 120; %mS/cm^2
cond_leak = 0.3; %mS/cm^2
Volt_K = -12; %mV
Volt_Na = 115; %mV
Volt_leak = 10.6; %mV
Volt_resting = 0; %mV %% set the resting potential at 0 and then re-adjust later to make it -70.
%when I just made the resting potential -70 the functions wouldn't work
%right. There would be action potentials when there were not supposed to be
%any
Cm = 1; %microFarads/cm^2

Vm = Volt_resting;
%initial conditions to variables that are changing with membrane
%voltage
alphan =0.01 * ( (10-Vm)/(exp((10- Vm)/10)-1) );
betan =.125 * exp(-Vm/80);
alpham = 0.1*( (25-Vm) / (exp( (25-Vm)/10) -1));
betam = 4*exp(-Vm/18);
alphah =0.07 * exp(-Vm/20);
betah = 1/(exp(( 30-Vm)/10)+1);
    
n = alphan/(alphan + betan);
m = alpham/(alpham + betam);
h = alphah/(alphah + betah);
    %no pulse in steady state
    I = 0;
    %this is needed to graph the conductance later on. I just had to make a vetor of
    %running conductances 
    gNa= m.^3*h*cond_Na;
    gK= n.^4*cond_K;
    
for t = 1:numel(time) - 1
    %t is an index of each step time for the whole time starting at 0 and
    %going to 10,000 (the total number of steps I am running the model for)
        
    %gating variables which change with changing Vm. 
    alphan(t) = 0.01 *( (10 - Vm(t))/(exp((10 - Vm(t))/10)-1));
    betan(t) = .125.*exp(-Vm(t)/80);
    alpham(t) = 0.1.*((25 - Vm(t))/(exp((25 -  Vm(t))/10)-1));
    betam(t) = 4.*exp(-Vm(t)/18);
    alphah(t) = 0.07 * exp(-Vm(t)/20);
    betah(t) =1/(exp((30 - Vm(t))/10)+1);

    %currents based on conductance and voltage
    INa = (m(t).^3).*cond_Na .* h(t).* (Vm(t) - Volt_Na);
    IK = (n(t).^4).*cond_K .* (Vm(t) - Volt_K);
    Ileak = cond_leak .* (Vm(t) - Volt_leak);
    Iion = I - IK - INa - Ileak; 

    % probability of Ion channels depending on gating variables that change
    % with a derivative function. Evaluate using Eulers Method
    m = [m (m(t)+step1*(alpham(t).*(1-m(t))-(betam(t)*m(t))))];
    n = [n (n(t)+step1*(alphan(t).*(1-n(t))-(betan(t)*n(t))))];
    h = [h  (h(t)+step1*(alphah(t).*(1-h(t))-(betah(t)*h(t))))];

    %re-evaluate the voltage potential and make it into a running vector to
    %plot
    Vm = [Vm Vm(t)+step1.*(Iion./Cm)]; 

    %make a vector of these conductances to graph later on
    gNa = [gNa (m(t).^3).*cond_Na*h(t)];
    gK = [gK (n(t).^4).*cond_K];
end
%re-adjust to make the membrane potential -70. I'm not sure why I couldn't
%make it -70 to begin with. The graph spikes and then settles to a little
%less than zero. This could have been due to the fact that a neuron cannot
%survive without some sort of action potential activity. Starting it at -70
%(I think) graphed what would happen if a neuron died due to a lack of
%activity. 
Vm = Vm - 70;
figure 
        plot(time, Vm)
        title('Membrane potential vs. time for steady state neuron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)') 
        
        figure
        c1 = plot(time,gNa);
        hold on
        c2 = plot(time, gK, 'r');
        legend([c1, c2], 'conductance for sodium', 'conductance for potassium')
        title('gNa and gK for steady state neuron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Coductance (mS/cm^2)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2  graphing with a step pulse
clear
%set the time that this function will run for
totaltime = 100; 
%set the step of each re-evaluation of the derivative functions
step1 = .01;
%make a vector of the times that this function will be re-evaluated (in
%milliseconds) starting at zero. I will use this later to plot my points
%and to run the loops.
time = 0:step1:totaltime;

%constants that were given
cond_K = 36; %mS/cm^2 
cond_Na = 120; %mS/cm^2
cond_leak = 0.3; %mS/cm^2
Volt_K = -12; %mV
Volt_Na = 115; %mV
Volt_leak = 10.6; %mV
Volt_resting = 0; %mV %% set the resting potential at 0 and then re-adjust later to make it -70.
%when I just made the resting potential -70 the functions wouldn't work
%right. There would be action potentials when there were not supposed to be
%any
Cm = 1; %microFarads/cm^2

Vm = Volt_resting;
%initial conditions to variables that are changing with membrane
%voltage
alphan =0.01 * ( (10-Vm)/(exp((10- Vm)/10)-1) );
betan =.125 * exp(-Vm/80);
alpham = 0.1*( (25-Vm) / (exp( (25-Vm)/10) -1));
betam = 4*exp(-Vm/18);
alphah =0.07 * exp(-Vm/20);
betah = 1/(exp(( 30-Vm)/10)+1);
    
n = alphan/(alphan + betan);
m = alpham/(alpham + betam);
h = alphah/(alphah + betah);
%no pulse in steady state
I = 0;
%this is needed to graph the conductance later on. I just had to make a vetor of
%running conductances 
gNa= m.^3*h*cond_Na;
gK= n.^4*cond_K;
    
for t = 1:numel(time) - 1
   %t is an index of each step time for the whole time starting at 0 and
    %going to 10,000 (the total number of steps I am running the model for)
        
    %gating variables which change with changing Vm. 
    alphan(t) = 0.01 *( (10 - Vm(t))/(exp((10 - Vm(t))/10)-1));
    betan(t) = .125.*exp(-Vm(t)/80);
    alpham(t) = 0.1.*((25 - Vm(t))/(exp((25 -  Vm(t))/10)-1));
    betam(t) = 4.*exp(-Vm(t)/18);
    alphah(t) = 0.07 * exp(-Vm(t)/20);
    betah(t) =1/(exp((30 - Vm(t))/10)+1);

    %currents based on conductance and voltage
    INa = (m(t).^3).*cond_Na .* h(t).* (Vm(t) - Volt_Na);
    IK = (n(t).^4).*cond_K .* (Vm(t) - Volt_K);
    Ileak = cond_leak .* (Vm(t) - Volt_leak);
        if t <50 %50 because that is the first 0.5 ms
        I = 5; %giving it 5microamps/cm^2 does not cause an action potential !!! that is so frustrating
        else 
        I = 0; %stop the pulse after 0.5 ms
        end
    
    Iion = I - IK - INa - Ileak; 

    % probability of Ion channels depending on gating variables that change
    % with a derivative function. Evaluate using Eulers Method
    m = [m (m(t)+step1*(alpham(t).*(1-m(t))-(betam(t)*m(t))))];
    n = [n (n(t)+step1*(alphan(t).*(1-n(t))-(betan(t)*n(t))))];
    h = [h  (h(t)+step1*(alphah(t).*(1-h(t))-(betah(t)*h(t))))];

    %re-evaluate the voltage potential and make it into a running vector to
    %plot
    Vm = [Vm Vm(t)+step1.*(Iion./Cm)]; 

    %make a vector of these conductances to graph later on
    gNa = [gNa (m(t).^3).*cond_Na*h(t)];
    gK = [gK (n(t).^4).*cond_K];
    

end
Vm = Vm - 70;
%re-adjust to make the membrane potential -70. I'm not sure why I couldn't
%make it -70 to begin with. The graph spikes and then settles to a little
%less than zero. This could have been due to the fact that a neuron cannot
%survive without some sort of action potential activity. Starting it at -70
%(I think) graphed what would happen if a neuron died due to a lack of
%activity. 
figure 
        plot(time, Vm)
        title('Membrane potential vs. time for pulse into a neuron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)') 
        
        figure
        c1 = plot(time,gNa);
        hold on
        c2 = plot(time, gK, 'r');
        legend([c1, c2], 'conductance for sodium', 'conductance for potassium')
        title('gNa and gK for a neuron given a pulse')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Coductance (mS/cm^2)')
        %% Part 3 graphing with a constant current of 5 microamps/cm^2
clear
%set the time that this function will run for
totaltime = 100; 
%set the step of each re-evaluation of the derivative functions
step1 = .01;
%make a vector of the times that this function will be re-evaluated (in
%milliseconds) starting at zero. I will use this later to plot my points
%and to run the loops.
time = 0:step1:totaltime;

%constants that were given
cond_K = 36; %mS/cm^2 
cond_Na = 120; %mS/cm^2
cond_leak = 0.3; %mS/cm^2
Volt_K = -12; %mV
Volt_Na = 115; %mV
Volt_leak = 10.6; %mV
Volt_resting = 0; %mV %% set the resting potential at 0 and then re-adjust later to make it -70.
%when I just made the resting potential -70 the functions wouldn't work
%right. There would be action potentials when there were not supposed to be
%any
Cm = 1; %microFarads/cm^2

Vm = Volt_resting;
%initial conditions to variables that are changing with membrane
%voltage
alphan =0.01 * ( (10-Vm)/(exp((10- Vm)/10)-1) );
betan =.125 * exp(-Vm/80);
alpham = 0.1*( (25-Vm) / (exp( (25-Vm)/10) -1));
betam = 4*exp(-Vm/18);
alphah =0.07 * exp(-Vm/20);
betah = 1/(exp(( 30-Vm)/10)+1);
    
n = alphan/(alphan + betan);
m = alpham/(alpham + betam);
h = alphah/(alphah + betah);
    %no pulse in steady state
    I = 0;
    %this is needed to graph the conductance later on. I just had to make a vetor of
    %running conductances 
    gNa= m.^3*h*cond_Na;
    gK= n.^4*cond_K;
    
    
for t = 1:numel(time) - 1
   %t is an index of each step time for the whole time starting at 0 and
    %going to 10,000 (the total number of steps I am running the model for)
        
    %gating variables which change with changing Vm. 
    alphan(t) = 0.01 *( (10 - Vm(t))/(exp((10 - Vm(t))/10)-1));
    betan(t) = .125.*exp(-Vm(t)/80);
    alpham(t) = 0.1.*((25 - Vm(t))/(exp((25 -  Vm(t))/10)-1));
    betam(t) = 4.*exp(-Vm(t)/18);
    alphah(t) = 0.07 * exp(-Vm(t)/20);
    betah(t) =1/(exp((30 - Vm(t))/10)+1);

    %currents based on conductance and voltage
    INa = (m(t).^3).*cond_Na .* h(t).* (Vm(t) - Volt_Na);
    IK = (n(t).^4).*cond_K .* (Vm(t) - Volt_K);
    Ileak = cond_leak .* (Vm(t) - Volt_leak);
    
    if t > 0 %for all times induce a current of 5 microamps into the neuron.
        I = 5; %5 microamps only causes one action potential !!!!!!! This was frustrating, a constant current should cause multiple action potentials
    else 
        I = 0;
    end
    Iion = I - IK - INa - Ileak; 

    % probability of Ion channels depending on gating variables that change
    % with a derivative function. Evaluate using Eulers Method
    m = [m (m(t)+step1*(alpham(t).*(1-m(t))-(betam(t)*m(t))))];
    n = [n (n(t)+step1*(alphan(t).*(1-n(t))-(betan(t)*n(t))))];
    h = [h  (h(t)+step1*(alphah(t).*(1-h(t))-(betah(t)*h(t))))];

    %re-evaluate the voltage potential and make it into a running vector to
    %plot
    Vm = [Vm Vm(t)+step1.*(Iion./Cm)]; 

    %make a vector of these conductances to graph later on
    gNa = [gNa (m(t).^3).*cond_Na*h(t)];
    gK = [gK (n(t).^4).*cond_K];
end
Vm = Vm - 70;
figure 
        plot(time, Vm)
        title('Membrane potential vs. time for constant 5 microamps/cm^2 into a neuron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)') 
        
        figure
        c1 = plot(time,gNa);
        hold on
        c2 = plot(time, gK, 'r');
        legend([c1, c2], 'conductance for sodium', 'conductance for potassium')
        title('gNa and gK for a neuron given a constant 5 microamps/cm^2 current')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Coductance (mS/cm^2)')
%%Part 4
%I want to show you an action potential caused by a constant large current of
%50microamps/cm^2 
        clear
%set the time that this function will run for
totaltime = 100; 
%set the step of each re-evaluation of the derivative functions
step1 = .01;
%make a vector of the times that this function will be re-evaluated (in
%milliseconds) starting at zero. I will use this later to plot my points
%and to run the loops.
time = 0:step1:totaltime;

%constants that were given
cond_K = 36; %mS/cm^2 
cond_Na = 120; %mS/cm^2
cond_leak = 0.3; %mS/cm^2
Volt_K = -12; %mV
Volt_Na = 115; %mV
Volt_leak = 10.6; %mV
Volt_resting = 0; %mV %% set the resting potential at 0 and then re-adjust later to make it -70.
%when I just made the resting potential -70 the functions wouldn't work
%right. There would be action potentials when there were not supposed to be
%any
Cm = 1; %microFarads/cm^2

Vm = Volt_resting;
%initial conditions to variables that are changing with membrane
%voltage
alphan =0.01 * ( (10-Vm)/(exp((10- Vm)/10)-1) );
betan =.125 * exp(-Vm/80);
alpham = 0.1*( (25-Vm) / (exp( (25-Vm)/10) -1));
betam = 4*exp(-Vm/18);
alphah =0.07 * exp(-Vm/20);
betah = 1/(exp(( 30-Vm)/10)+1);
    
n = alphan/(alphan + betan);
m = alpham/(alpham + betam);
h = alphah/(alphah + betah);
    %no pulse in steady state
    I = 0;
    %this is needed to graph the conductance later on. I just had to make a vetor of
    %running conductances 
    gNa= m.^3*h*cond_Na;
    gK= n.^4*cond_K;
    
    
for t = 1:numel(time) - 1
   %t is an index of each step time for the whole time starting at 0 and
    %going to 10,000 (the total number of steps I am running the model for)
        
    %gating variables which change with changing Vm. 
    alphan(t) = 0.01 *( (10 - Vm(t))/(exp((10 - Vm(t))/10)-1));
    betan(t) = .125.*exp(-Vm(t)/80);
    alpham(t) = 0.1.*((25 - Vm(t))/(exp((25 -  Vm(t))/10)-1));
    betam(t) = 4.*exp(-Vm(t)/18);
    alphah(t) = 0.07 * exp(-Vm(t)/20);
    betah(t) =1/(exp((30 - Vm(t))/10)+1);

    %currents based on conductance and voltage
    INa = (m(t).^3).*cond_Na .* h(t).* (Vm(t) - Volt_Na);
    IK = (n(t).^4).*cond_K .* (Vm(t) - Volt_K);
    Ileak = cond_leak .* (Vm(t) - Volt_leak);
    
    if t > 0 %for all times induce a current of 5 microamps into the neuron.
        I = 50; %5 microamps only causes one action potential !!!!!!! This was frustrating, a constant current should cause multiple action potentials
    else 
        I = 0;
    end
    Iion = I - IK - INa - Ileak; 

    % probability of Ion channels depending on gating variables that change
    % with a derivative function. Evaluate using Eulers Method
    m = [m (m(t)+step1*(alpham(t).*(1-m(t))-(betam(t)*m(t))))];
    n = [n (n(t)+step1*(alphan(t).*(1-n(t))-(betan(t)*n(t))))];
    h = [h  (h(t)+step1*(alphah(t).*(1-h(t))-(betah(t)*h(t))))];

    %re-evaluate the voltage potential and make it into a running vector to
    %plot
    Vm = [Vm Vm(t)+step1.*(Iion./Cm)]; 

    %make a vector of these conductances to graph later on
    gNa = [gNa (m(t).^3).*cond_Na*h(t)];
    gK = [gK (n(t).^4).*cond_K];
end
Vm = Vm - 70;
figure 
        plot(time, Vm)
        title('Membrane potential vs. time for constant 50 microamps/cm^2 into a neuron')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)') 
        
        figure
        c1 = plot(time,gNa);
        hold on
        c2 = plot(time, gK, 'r');
        legend([c1, c2], 'conductance for sodium', 'conductance for potassium')
        title('gNa and gK for a neuron given a constant 50 microamps/cm^2 current')
        axis([0,100,-100,100])
        xlabel('Time (ms)')
        ylabel('Coductance (mS/cm^2)')
        
        %The action potentials are beautiful. It took so long to get these
        %that if I only showed you the 5microamp current it wouldn't have
        %showed everything I modeled. The first action potential is larger
        %than the later ones, which is what we learned in class. There is a
        %clear refractory period. On the conductance graph you can clearly
        %see the sodium spiking before the potassium, which shows the
        %sodium going into the cell, and then the potassium starting slower
        %and going out of the cell. Beautiful. 
