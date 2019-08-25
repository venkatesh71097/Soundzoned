[x1, Fs1] = audioread('ica1.wav');
[x2, Fs2] = audioread('ica2.wav');

m = size(x1,1); % size of each signal
n = 2; % Number of sound sources
A = randn(n, n); %  Random 2 X 2 mixing matrix
x = A*[x1';x2']; % Mixed signals
c = cov(x');
sq = inv(sqrtm(c));
mx = mean(x, 2)';
x = x - mx'*ones(1, size(x, 2));
x = sq*x;
w1 = randn(n, 1); % Random weight vector
w1 = w1/norm(w1,2); % make 'w' a unit vector
w0 = randn(n, 1);
w0 = w0/norm(w0, 2);
while abs(abs(w0'*w1)-1) > 0.000001
    w0 = w1;
    w1 = x*G(w1'*x)'/m - sum(DG(w1'*x))*w1/m; % This step is supposed to perform:
                                          % w = E{xg(w^{T}*x)} - E{g'(w^{T}*x)}w
    w1 = w1/norm(w1, 2);
end

w2 = randn(n, 1);
w2 = w2/norm(w2,2);
w0 = randn(n, 1);
w0 = w0/norm(w0, 2);

while abs(abs(w0'*w2)-1) > 0.001
    w0 = w2;
    w2 = x*G(w2'*x)'/m - mean(DG(w2'*x), 2)*w2;
    w2 = w2 - w2'*w1*w1;
    w2 = w2/norm(w2, 2);
end
w = [w1 w2];
s = w*x;
s1 = s(1,:);
s2 = s(2,:);
plot(x1);
figure;
plot(x2);
figure;
plot(s1);
figure;
plot(s2);

audiowrite('output1.wav',s1,44100);
audiowrite('output2.wav',s2,44100); 

