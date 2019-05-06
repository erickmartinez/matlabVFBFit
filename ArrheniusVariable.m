function [x] = ArrheniusVariable(x0,Ea,T)
    %     k = 1.38064852; % E-23 J/K
    %     eV = 1.60217662; % E-19 J
    kT = T*(1.38064852/1.60217662)*1E-4;
    x = x0*exp(-Ea/kT);
end