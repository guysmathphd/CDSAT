function [ff,P1,P1a] = myfft(all_times_sorted,total_weight)
all_times_sorted_seconds = seconds(all_times_sorted - all_times_sorted(1));
[m,i] = max(diff(all_times_sorted_seconds));
T = m;
Fs = 1/T;
t = 0:T:all_times_sorted_seconds(end);
L = length(t);
if mod(L,2) ~= 0
%     all_times_sorted_seconds = all_times_sorted_seconds(1:end-1);
    L = L - 1;
    t = t(1:end-1);
end

ind1 = 1;ind2 = 1;
while ind1 <= length(t)
    if all_times_sorted_seconds(ind2) == t(ind1)
        total_weight_sampled(ind1) = total_weight(ind2);
        ind1 = ind1+1; ind2 = ind2+1;
    elseif all_times_sorted_seconds(ind2) < t(ind1)
        ind2 = ind2+1;
    elseif all_times_sorted_seconds(ind2) > t(ind1)
        total_weight_sampled(ind1) = total_weight(ind2-1);
        ind1 = ind1 + 1;
    end
end
Y = fft(total_weight_sampled);
P2 = abs(Y/L);P2a = angle(Y);
P1 = P2(1:L/2+1);P1a = P2a(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
ff = Fs*(0:(L/2))/L;

end