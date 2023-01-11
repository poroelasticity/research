function N_interp = Ny()

        % a set of discrete data points y ranging from y=0 to y=2e8
        number_data_points =  500;
        y_end              =  2e8;
        % the sampling points are allocated in a way that the key trends in \mathrm{N}\left(y\right) can be well captured.
        y                  =  [0,exp(log(0.001):(log(y_end)-log(0.001))/(number_data_points-2):log(y_end))];       
        N                  =  y;

        % when 0<y<=2000, results of N(y) are calculated by numerically 
        % integrating the integrals in M0(y) and M1(y);
        y_1 = y((y<=2000));

        for i=1:length(y_1)

              ft = @(t) (1-t.^2).^(1/2).*exp(-y_1(i)*t);
              M1 = -2/pi*y_1(i)*integral(ft,0,1);

              ft = @(t) exp(-y_1(i)*cos(t));
              M0 = -2/pi*integral(ft,0,pi/2);

            N(i) = 1/y_1(i)*(M0 - 2*M1/y_1(i)); 
            
        end

        % when y=0, N(0) equal to 4/3/pi

        N(1) = 2/3/pi;
        
        % when y>2000, N(y) = 2*y^(-2)/pi
        N((y>2000)) = 2/pi*(y(y>2000)).^(-2);


        % create a function that linearly interpolates between vectors y 
        % and N(y) using the MATLAB command interp1
        N_interp = @(y_N) interp1(y,N,y_N);

end