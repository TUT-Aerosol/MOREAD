function [out] = CPCsigmoid(Dp,cutoff)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

out = 1./sqrt(1+exp(-6e9.*(Dp-cutoff)));

end

