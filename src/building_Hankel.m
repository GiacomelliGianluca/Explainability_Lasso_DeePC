function H = building_Hankel(Data, T_h, d, n_x, flag_PE)
%Function that constructs Hankel matrix

% T_h: prediction horizon. It is composed by the past horizon T_ini and the
% future horizon L (T_h = T_ini + L). The initial instant assumed for the
% Hankel matrix is t = 0
% d: data dimensionality (corresponing to the number of elements for an istant t)
% n_x: state dimensionality
% flag_PE: 1: test the persistency of excitation 0: no test


N_Data = length(Data);

% Persistency of excitation (PE) test
if flag_PE
    PE_check(Data, d*(T_h + n_x));
end




%Partitioning raws and columns
raw_data = Data(1:d*T_h);
column_data = Data(d*T_h:N_Data);

H = hankel(raw_data, column_data);



end