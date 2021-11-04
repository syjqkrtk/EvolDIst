seqlen = zeros(1,25);

for i = 1:25
    seqlen(i) = length(NewSequence3{i});
end

plot(0:0.01:0.25,length(sequence)*ones(1,26));
hold on
plot(0:0.01:0.25,[length(sequence) seqlen]);
xlabel('t');
ylabel('sequence length');
legend('Raw sequence', 'Simulated sequence');
ylim([0 50000]);