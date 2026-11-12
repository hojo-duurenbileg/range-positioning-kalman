function [ Qk ] = QkCalc( deltaT, Sp, Sv )


Qk = [(Sp * deltaT + (1/3)*Sv*deltaT^3) 0 (0.5*Sv*deltaT^2) 0;
    0 (Sp * deltaT + (1/3)*Sv*deltaT^3) 0 (0.5*Sv*deltaT^2);
    (0.5*Sv*deltaT^2) 0 Sv*deltaT 0;
    0 (0.5*Sv*deltaT^2) 0 Sv*deltaT];


end

