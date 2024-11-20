function I = trapeze(signal, dt, cdtinit)
    % Trapezoidal integration of a signal
    %
    % Inputs:
    %   signal - input signal to be integrated
    %   dt - time step
    %   cdtinit - initial condition
    %
    % Output:
    %   I - integrated signal

    n = length(signal);
    I = zeros(1, n);
    I(1) = cdtinit;

    for i = 2:n
        I(i) = I(i-1) + (dt/2) * (signal(i) + signal(i-1));
    end

  end