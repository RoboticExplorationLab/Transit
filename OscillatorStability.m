samp_rate = 500e3;
t = (0:(samp_rate-1))/samp_rate;
num_steps = floor(length(raw_data)/samp_rate);

fhist = zeros(num_steps,1);
for k = 1:num_steps
    inds = (1:samp_rate)+(k-1)*samp_rate;
    samp = raw_data(inds);
    [freq,fft] = easyfft(t,samp,samp_rate/5,1);
    [~,maxind] = max(fft(3:end));
    fhist(k) = freq(maxind);
end

figure()
plot(fhist-mean(fhist));

xlabel('Time (sec)');
ylabel('Frequency Drift (Hz)');