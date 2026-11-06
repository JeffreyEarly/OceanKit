N=[50; 100; 200; 400; 600; 800; 1000; 1200];
solverTime = zeros(size(N));

for i=1:length(N)
    A = rand(N(i));
    B = rand(N(i));
    tic
    [~,~] = eig( A, B );
    solverTime(i) = toc;
end

%%
figure
plot(N,solverTime,LineWidth=2), ylog, xlog
p = polyfit(log(N),log(solverTime),1);
model = @(N) exp(p(2))*N.^p(1);
hold on, plot(N,model(N),Color=0*[1 1 1],LineWidth=1)

%% What is the optimal search increment given how long it takes to compute an EVP?

startingGuess = 1400;
trueAnswer = 1500;

bestIncrement = 1;
bestCumulativeTime = Inf;

for dEVP = 1:startingGuess
    cumulativeTime = 0;
    currentGuess = startingGuess;
    while currentGuess < trueAnswer+2*dEVP % you won't know you crossed the true answer until you got there AND exceeded it with your convergence test
        cumulativeTime = cumulativeTime + model(currentGuess);
        currentGuess = currentGuess + dEVP;
    end

    if cumulativeTime < bestCumulativeTime
        bestCumulativeTime = cumulativeTime;
        bestIncrement = dEVP;
    end
end
disp("best cumulative time: " + string(bestCumulativeTime) + " with increment " + bestIncrement)

%% Test other increments and see how much longer they take
cumulativeTime = 0;
dEVP = 100;
currentGuess = startingGuess;
while currentGuess < trueAnswer+2*dEVP % you won't know you crossed the true answer until you got there AND exceeded it with your convergence test
    cumulativeTime = cumulativeTime + model(currentGuess);
    currentGuess = currentGuess + dEVP;
end
disp("cumulative time: " + string(cumulativeTime))

% Roughly, it seems that you might want to just start really low and take
% big-ish increments?