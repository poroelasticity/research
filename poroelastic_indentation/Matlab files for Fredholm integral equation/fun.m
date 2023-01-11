function results = fun(Gauss_x,x12_M_interp1,x12_M_interp2,interval)
    
    Gauss_x1 = Gauss_x(find(Gauss_x <= 1/interval));
    Gauss_x2 = Gauss_x(find(Gauss_x >= 1/interval));

    results = [x12_M_interp1(Gauss_x1),x12_M_interp2(Gauss_x2)];

end

