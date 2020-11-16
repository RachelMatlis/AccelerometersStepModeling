function AccelStepModelingAlgo1
    close all;
    format long
    % Load h5 file in MATLAB
    AccelSignal = h5read('20140126-134055-Base.h5','/SI-001007/Calibrated/Accelerometers');
   %RA = h5read('20140126-134055-Base.h5','/SI-001011/Calibrated/Accelerometers');

    Fs = 128;
    Ts = 1/Fs;
    
    axisY = AccelSignal(1,:);
    axisX = AccelSignal(3,:);
    axisZ = AccelSignal(2,:);
    
    N = length(axisY);
    t =linspace(0, N * Ts, N);
    
    %Subtracting the mean from the data will remove any constant effects, such as gravity
    axisY = axisY - mean(axisY);
    axisX = axisX - mean(axisX);
    axisZ = axisZ - mean(axisZ); 
        
    AccelSignal(1,:) = axisY;
    AccelSignal(2,:) = axisZ;
    AccelSignal(3,:) = axisX;

% Divide all of the accelerometer signals to separate steps using the height accelerometer

    %Apply LPF to to smooth out the signal     
    smoothedY = lpfFilter(axisY, 4550);
    
   %Each peak in acceleration corresponds to a step being taken while walking
    periodTime = getPeriodTime(smoothedY, Fs);
    minPeakHeight = std(smoothedY);
    minPeakDistance = Fs * periodTime;
    [pks, locs] = findpeaks(smoothedY, 'MinPeakDistance',  minPeakDistance, 'MinPeakHeight', minPeakHeight);
   
    %Detect step start and end:
    for i = 1 : length(locs)
           locs(i) = getZeroPointBeforeLocalMaxima(minPeakDistance, locs(i), smoothedY);
    end
    
    %Build steps matrix
    steps = BuildStepsVectorMatrix(locs, minPeakDistance);
    
    %Implement interpolation in order to make each step signal samples equal.
    interpolatedSteps = stepsInterpolation(steps, AccelSignal);
    
% Calculate the mean step and its variance.
    maxStepLength = ceil(max(steps(:, 2) - steps(:, 1)));
    meanStep = zeros(3, maxStepLength);
    meanStep(1, :) = mean(interpolatedSteps(:,1,:)); 
    meanStep(2, :) = mean(interpolatedSteps(:,2,:)); 
    meanStep(3, :) = mean(interpolatedSteps(:,3,:));
    
% Calculate var and std for each axis mean step.
    varY = permute(var(interpolatedSteps(:,1,:)), [3 2 1]);
    varZ = permute(var(interpolatedSteps(:,2,:)), [3 2 1]);
    varX = permute(var(interpolatedSteps(:,3,:)), [3 2 1]);
    
    stdY = permute(std(interpolatedSteps(:,1,:)), [3 2 1]); 
    stdZ = permute(std(interpolatedSteps(:,2,:)), [3 2 1]);
    stdX = permute(std(interpolatedSteps(:,3,:)), [3 2 1]);
    
% Project Plots:

    %Plottig the signals with the marks of the start of each step
    
    %axis Y
    figure()
    subplot(3,1,1)
    plot(t, axisY);
    hold on
    plot(t,smoothedY,'r');
    plot(t(locs),axisY(locs), 'k', 'Marker', 'o', 'LineStyle', 'none');
    title('Axis Y');
    legend('Signal','Low Pass Filter', 'Step Start');
    hold off
    grid on
    
    %axis Z
    subplot(3,1,2)
    plot(t,axisZ);
    hold on
    plot(t(locs), axisZ(locs), 'r', 'Marker', 'o', 'LineStyle', 'none');
    title('Axis Z');
    hold off
    grid on
    
    %axis X
    subplot(3,1,3)
    plot(t,axisX);
    hold on
    plot(t(locs), axisX(locs), 'r', 'Marker', 'o', 'LineStyle', 'none');
    title('Axis X');
    hold off
    grid on
    
% Plot the mean step with its std for all of the accelerometer.

    %Plot all the steps axis in blue and the mean step in red
    
    figure();
    subplot(3,1,1);
    hold on
    for i = 1 : length(steps)
        plot(permute(interpolatedSteps(i,1,:), [3 2 1]), 'blue');
    end
    plot(meanStep(1,:),'red', 'LineWidth',2);
    title('Axis Y- After interpolation and the Mean Step');
    hold off

    subplot(3,1,2);
    hold on
    for i = 1 : length(steps)
        plot(permute(interpolatedSteps(i,2,:), [3 2 1]), 'blue');
    end
    plot(meanStep(2,:),'red', 'LineWidth',2);
    title('Axis Z- After interpolation and the Mean Step');
    hold off
    
    subplot(3,1,3);
    hold on
    for i = 1 : length(steps)
        plot(permute(interpolatedSteps(i,3,:), [3 2 1]), 'blue');
    end
    plot(meanStep(3,:),'red', 'LineWidth',2);  
    title('Axis X- After interpolation and the Mean Step');
    hold off
        
% Plot the Mean Step with it's Variance and Standard deviation
   
    figure();
    plot(meanStep(1,:));
    hold on
    plot(varY, '--');
    plot(stdY, '--');
    legend('Axis Y signal', 'Variance', 'Standard Deviation');
    hold off
    
    figure();
    plot(meanStep(2,:));
    hold on
    plot(varZ, '--');
    plot(stdZ, '--');
    legend('Axis Z signal', 'Variance', 'Standard Deviation');
    hold off
    
    figure();
    plot(meanStep(3,:));
    hold on
    plot(varX, '--');
    plot(stdX, '--');
    legend('Axis X signal', 'Variance', 'Standard Deviation');
    hold off
    
end

function Signal= lpfFilter(signal, cutOfInd)

    editedSignal = fftshift(fft(signal));

    windowSize = length(signal);
    cutOfRange = zeros(size(signal));
    lowCenter = floor(windowSize/2);
    highCenter = ceil(windowSize/2);

    cutOfRange(lowCenter - cutOfInd :  highCenter + cutOfInd) = 1;
    lpf = editedSignal.* cutOfRange;
    
    Signal = real(ifft(ifftshift(lpf))); 
end

function periodicityTime = getPeriodTime(signal, fs)
    signalNorm = signal-mean(signal);
    [autocor,lags] = xcorr(signalNorm,ceil(fs * (length(signal) / fs / 10)),'coeff');
    [pksh,lcsh] = findpeaks(autocor, 'MinPeakHeight', 0.05);
    periodicityTime = mean(diff(lcsh))/fs;
end

function groundZeroPoint = getZeroPointBeforeLocalMaxima(periodValuesCount, location, signal)

    moveBackLocation = location - periodValuesCount;
    while (abs(signal(location)) >= 0 && location >= moveBackLocation)
       location = location - 1;
    end
    
    groundZeroPoint = location;   
end

function steps = BuildStepsVectorMatrix(stepsLocationVec, minPeakDis)

    rows = length(stepsLocationVec);
    column = 2;
    steps = zeros(rows,column);
    for i = 1 : rows - 1
        if stepsLocationVec(i+1) > (stepsLocationVec(i) + minPeakDis * 1.2)
            steps(i,:) = [stepsLocationVec(i), stepsLocationVec(i) + minPeakDis];
        else
            steps(i,:) = [stepsLocationVec(i), stepsLocationVec(i+1) - 1];
        end
    end
    
    lastStep = length(steps);
    steps(lastStep,:) = [stepsLocationVec(lastStep), stepsLocationVec(lastStep) + minPeakDis];
end

function interpolatedSteps = stepsInterpolation(steps, AccelSignal)

    maxStepLength = max(steps(:, 2) - steps(:, 1));
    stepsSignal = zeros(length(steps), 3, ceil(maxStepLength));
    
    for i = 1 : length(steps)
        curStepIndex = 0:(steps(i,2)-steps(i,1)); 
        curStepValueY = AccelSignal(1, steps(i,1):steps(i,2));
        curStepValueZ = AccelSignal(2, steps(i,1):steps(i,2)); 
        curStepValueX = AccelSignal(3, steps(i,1):steps(i,2));

        interpStepIndex = 0 : (length(curStepIndex)/maxStepLength):  length(curStepIndex) - 0.1; 
        stepsSignal(i, 1,:) = interp1(curStepIndex, curStepValueY, interpStepIndex, 'PCHIP'); 
        stepsSignal(i, 2,:) = interp1(curStepIndex, curStepValueZ, interpStepIndex, 'PCHIP'); 
        stepsSignal(i, 3,:) = interp1(curStepIndex, curStepValueX, interpStepIndex, 'PCHIP'); 
    end 
    
    interpolatedSteps = stepsSignal;
end
