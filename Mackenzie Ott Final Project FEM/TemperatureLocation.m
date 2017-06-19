function Temp = TemperatureLocation(T, s,t)
Temp = (1/4)*(T(3) *(1+s)*(1+t) + T(1) *(1-s)*(1+t) + T(2)*(1-t)*(1-s) + T(4)*(1+s)*(1-t));