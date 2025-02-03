function [sumKahan, compensation] = kahanSumIterative(sumKahan, compensation, number)
    y = number - compensation;
    t = sumKahan + y;
    compensation = (t - sumKahan) - y;
    sumKahan = t;
end
