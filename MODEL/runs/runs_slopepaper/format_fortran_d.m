function[out] = format_fortran_d(in)




if length(in)>1,
    if isvector(in)
        out = '';
        for i = 1:length(in)
            if(in(i)~=0)
                b1 = in(i)./(10.^(floor(log10(in(i)))));
                b2 = floor(log10(in(i)));
                if b2>=0,
                    out=[out sprintf('%6.5fE+%i ',b1,b2)]; % strcat does not work because of trailing space
                else
                    out=[out sprintf('%6.5fE%i ',b1,b2)];
                end
            else
                out=[out sprintf('%6.5fE+%i ',0,0)];
            end
        end
    else
        error('matrix values cannot be formatted like this!')
    end
    
else
    if (in~=0)
        b1 = in./(10.^(floor(log10(in))));
        b2 = floor(log10(in));
        if b2>=0,
            out=sprintf('%6.5fE+%i ',b1,b2);
        else
            out=sprintf('%6.5fE%i ',b1,b2);
        end
    else
        out=sprintf('%6.5fE+%i ',0,0);
    end
end
end