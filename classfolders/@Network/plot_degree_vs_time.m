function obj = plot_degree_vs_time(obj)
%net05a
total_weight = General.load_var(fullfile(obj.path,'total_weight_vs_all_times_sorted'));
all_times_sorted = General.load_var(fullfile(obj.path,'all_times_sorted'));

name = 'net05a';figdesc = 'Total network weight vs time';
f = figure('Name',name,'NumberTitle','off');
plot(all_times_sorted,total_weight);
xlabel('Times');
ylabel('$k_{total} = \sum_i k_{i}$','Interpreter','latex');
title({[name ' ' obj.name ' ' obj.desc];figdesc});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));

%net05b
all_times_sorted_seconds = seconds(all_times_sorted - all_times_sorted(1));
L = all_times_sorted_seconds(end);
if mod(L,2) ~= 0
%     all_times_sorted_seconds = all_times_sorted_seconds(1:end-1);
    L = L - 1;
end
Fs = 1; T = 1;
t = (0:L-1)*T;total_weight_sampled = zeros(size(t));ind = 1;
i = 0;
while i <= t(end)
%     ind = find(all_times_sorted_seconds >= i,1,'first');
    if all_times_sorted_seconds(ind) == i
        total_weight_sampled(i+1) = total_weight(ind);
        i = i + 1; ind = ind + 1;
    elseif all_times_sorted_seconds(ind) < i
        ind = ind+1;
    elseif all_times_sorted_seconds(ind) > i
        total_weight_sampled(i+1) = total_weight(ind-1);
        i = i + 1;
    end
end
Y = fft(total_weight_sampled);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L*24*3600;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
end