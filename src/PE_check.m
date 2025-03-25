function flag_PE = PE_check(Data,depth)
%PE_CHECK checks the persistency of excitation of a signal


H_PE = hankel(Data(1:depth), Data(depth:end));
r_H_PE = rank(H_PE);

if r_H_PE ~= (depth)
    keyboard
    flag_PE = 0;
else
    flag_PE = 1;
end
end

