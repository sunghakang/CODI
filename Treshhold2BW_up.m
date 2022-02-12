function u = Treshhold2BW_up (x,tresh)
[n,m]=size(x);
for i=1:n
    for j=1:m
        if x(i,j)<tresh
            u(i,j)=0;
        else
            u(i,j)=255;
        end
    end
end
end