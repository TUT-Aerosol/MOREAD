function[out] = format_fortran_d(in)

b1 = in./(10.^(floor(log10(in))));
b2 = floor(log10(in));
if b2>=0,
    out=sprintf('%6.5fE+%i ',b1,b2);
else
    out=sprintf('%6.5fE%i ',b1,b2);
end
end