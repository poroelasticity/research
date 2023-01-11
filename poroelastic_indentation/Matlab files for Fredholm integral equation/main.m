clear;clc;

% this code is for case I surface drainage condition

%% Step 1: construct an interpolation function to evaluate N(s,x,m)
N_interp = Ny() ;

%% step 2: evaluate x12_Msx through numerical integration
% create a set of sampling data points from 0 to inf
interval      =  0.0025                   ; % segment length
x_0_1         =  0:interval:1             ; % In [0,1], the points are evenly spaced, with the length of each segment being 0.0025.
x_1_interval  =  1:-interval:interval     ;
x_1_inf       =  1./x_1_interval          ; % In [1,inf), the points are distributed by mapping the points in [0,1] to [1,inf) via 1/x. 
x_0_inf       =  [x_0_1(1:end-1),x_1_inf] ;  

% create sufficient number of Gauss points in [0,1]. 
% the subroutine file named "lgwt.m” is used to calculate the location of Gauss points and their weights. 
Gauss_m = [] ;
Gauss_w = [] ;
for i = 1:length(x_0_1)-1
    [Gauss_m_sub,Gauss_w_sub]  =  lgwt(32,x_0_1(i),x_0_1(i+1)) ;
    Gauss_m = [Gauss_m, Gauss_m_sub']                          ; 
    Gauss_w = [Gauss_w, Gauss_w_sub']                          ; 
end

% determine the value of x12_Msx at each discrete point through a weighted sum.
s = 1; % Laplace coefficient
for i = 1:length(x_0_inf)
    x          = x_0_inf(i); 
    x12_N      = sqrt(s)*(N_interp(sqrt(s)*abs(x-Gauss_m))+N_interp(sqrt(s)*(x+Gauss_m)));
    % Computation of x12_Nsxm is then based on the interpolation function, which is more efficient because it avoids numerical integration.
    x12_M(i)   = sum((1-Gauss_m.^2).*x12_N.*Gauss_w); % paraboloidal: 1-Gauss_m.^2, conical: 1-Gauss_m, cylindrical: 1
end
% the following functions will be used to determine a1 in the next step
x12_M_interp1 = @(x) interp1(x_0_inf,x12_M,x); % create a function that interpolates between vectors x and x12_Msxd
coeff         = x12_M(end)/x_0_inf(end)^(-2);
x12_M_interp2 = @(x) coeff*x.^(-2);            % create a function to determine x12_Msx beyond x = 400 by using its asymptotic expressions


%% step 3: evaluate x12_asx through numerical integration
%create sufficient number of Gauss points in [1,800]. 
m_1_inf_aug = [x_1_inf,2/interval]; % replace the upper bound with a large value of 800
Gauss_m   = [];
Gauss_w   = [];
Gauss_num = 256*2;
for i = 1:length(m_1_inf_aug)-1
    [Gauss_m_sub,Gauss_w_sub]  =  lgwt(Gauss_num,m_1_inf_aug(i),m_1_inf_aug(i+1));
    Gauss_m = [Gauss_m, Gauss_m_sub']; 
    Gauss_w = [Gauss_w, Gauss_w_sub']; 
end

x12_a_storage = zeros(19,length(x_0_inf)); % values of all a_n (n>=1) will be stored in a matrix

% determine the value of x12_asx at each discrete point by a weighted sum.

for i_omega = 1:19

    i_omega

    for i = 1:length(x_0_inf)

        x          =  x_0_inf(i); 

        x12_N   =  sqrt(s)*(N_interp(sqrt(s)*abs(x-Gauss_m))+N_interp(sqrt(s)*(x+Gauss_m)));

        if i_omega == 1

            x12_M_1_inf_gauss  = fun(Gauss_m,x12_M_interp1,x12_M_interp2,interval); % determines the value of x12_M
            
            x12_a(i)  = sum(x12_N.*x12_M_1_inf_gauss.*Gauss_w); % for a1
        
        else
            
            x12_a_1_inf_gauss  = fun(Gauss_m,x12_a_interp1,x12_a_interp2,interval); % determines the value of x12_a

            x12_a(i)  = sum(x12_N.*x12_a_1_inf_gauss.*Gauss_w); % for a2, a3, a4, a5 ...
        
        end    
    end
    
    x12_a_storage(i_omega,:) = x12_a; % values of all a_n (n>=1) is stored in a matrix

    % the following functions will be used to determine an+1 in the next step
    x12_a_interp1 = @(x) interp1(x_0_inf,x12_a,x); % create a function that interpolates between vectors x and x12_asx
    coeff         = x12_a(end)/x_0_inf(end)^(-2);
    x12_a_interp2 = @(x) coeff*x.^(-2);            % create a function to determine x12_asx beyond x = 400 by using its asymptotic expressions
    
    
end 

%% step 4: evaluate x12_theta1
% take the sum of all a_n, when omega = or < 0.5, the series converges rather fast
% the partial sum with only 20 terms is sufficient to give a satisfactory approximation 
omega = 0.5;
x12_theta1 = omega*x12_M;
for i_omega = 1:19
    x12_theta1 = x12_theta1 - (-omega)^(i_omega + 1)*x12_a_storage(i_omega,:);
end

