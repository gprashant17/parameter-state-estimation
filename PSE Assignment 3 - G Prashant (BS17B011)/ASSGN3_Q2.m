%% Question 2a
N = 1000;  %Observations
R = 300;   %Realizations
n = 11+1;  %Roll no: BS17B011 -> n = 11 + 1
fn = n/N;  %Frequency

v_mat = randn(R,N); %GWN matrix
cos_vector = cos(2*pi*fn*0:(N-1));
sin_vector = sin(2*pi*fn*0:(N-1));

%Fourier Coefficients an and bn
an = sum(v_mat.*cos_vector,2);
bn = sum(-v_mat.*sin_vector,2);

%Plotting Histogram
figure;
histogram(an)
xlabel("an")
ylabel("f(an)")
title("Histogram of an")
figure;
histogram(bn)
xlabel("bn")
ylabel("f(bn)")
title("Histogram of bn")

%%KS Test for Normality Test
x = an/sqrt(sum(cos_vector.^2));
[h,p] = kstest(x);
if h==0
    fprintf("an is Gaussian Distributed\n")
else
    fprintf("an is not Gaussian Distributed\n")
end

y = bn/sqrt(sum(sin_vector.^2));
[h,p] = kstest(y);
if h==0
    fprintf("bn is Gaussian Distributed\n")
else
    fprintf("bn is not Gaussian Distributed\n")
end
Pn = (an.^2 + bn.^2)/N;
%% Question 2d MC simluation for Consistency of P(fn)
conv = zeros(1,length(1:500:1e5));
i = 1;
for N = 1:500:1e5
    v_mat = randn(R,N);
    cos_vector = cos(2*pi*fn*0:(N-1));
    sin_vector = sin(2*pi*fn*0:(N-1));
    an = sum(v_mat.*cos_vector,2);
    bn = sum(-v_mat.*sin_vector,2);
    Pn = (an.^2 + bn.^2)/N;
    conv(i) = var(Pn);
    i = i + 1;
end

plot(1:500:1e5,conv,"-")
xlabel("N")
ylabel("MSE")
title("Mean Square sense consistency of P(f)")