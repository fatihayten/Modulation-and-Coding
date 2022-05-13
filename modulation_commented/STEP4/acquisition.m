function [n, df] = acquisition(received_signal, pilot_symbols, T_symbol, K_parameter)
%{
    input
    rx: received symbols
    p:  pilot symbols
    T:  symbol duration
    K:  K-factor

    output
    n:  arrival time estimate
    df: cfo estimate
%}

N = length(pilot_symbols);
L = length(received_signal);
D = zeros(K_parameter, L-N+1);

for k = 1 : K_parameter
    for n = 1 : L - N + 1
        tmp = 0;
        for l = k : N-1
            tmp1 = conj(received_signal(n+l)) * pilot_symbols(l+1);
            
            tmp2 = conj(received_signal(n+l-k)) * pilot_symbols(l-k+1);
            
            tmp = tmp + tmp1 * conj(tmp2);
        end
        D(k,n) = tmp;
        
    end 
    D(k,:) = D(k,:)/(N-k);
end

% Time of arrival estimate
tmp = sum(abs(D),1);
[~, n] = max(tmp);

% CFO estimate
k = (1:K_parameter).' ;
df = sum(angle(D(k,n))./(2*pi*k*T_symbol), 1);
df = - df/K_parameter;

end