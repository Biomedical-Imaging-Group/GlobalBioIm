function len = fft_best_dim(len)
%%  fft_best_dim(len);
%     Return the smallest integer which is greater or equal LEN and which is
%     a multiple of powers of 2, 3 and/or 5

best= 2*len;
i5 = 1;
while true
    i3 = i5;
    while true
        i2 = i3;
        while i2 <len
            i2= i2*2;
        end
        if i2 ==len
            return;
        end
        if (i2-len < best-len)
            best= i2;
        end
        i3 = i3*3;
        if i3>len
            break;
        end
    end
    
    i5=i5*5;
    if i5>len
        break;
    end
end
len = best;
end