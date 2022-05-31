function F = in_transducer(M,a)
%     a = 1;
    t1 = 22;
    t0 = t1;
    t2 = 0.5;
    %t2 = 0.0335;
    F = zeros(size(M));
    [m,n] = size(M);
    for i = 1:m
        for j = 1:n
            if M(i,j) > 0
                if M(i,j) > t1
                    F(i,j) = t2*(M(i,j)/t1)^2;
                else
                    F(i,j) = t2*(M(i,j)/t1)^(1/a);
                end
            else
                M(i,j) = abs(M(i,j));
                if M(i,j) > t1
                    F(i,j) = -1*t2*(M(i,j)/t0)^2;
                else
                    F(i,j) = -1*t2*(M(i,j)/t0)^(1/a);
                end
            end
        end
    end
end