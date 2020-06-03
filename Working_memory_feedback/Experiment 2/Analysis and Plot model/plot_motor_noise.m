% Plot the motor noise from control experiment
subjectID = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh', 'average'}; % 'll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh', 
motor_noise = NaN(1, length(subjectID));    
for ii = 1 : length(subjectID)
    [binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForFitting(subjectID(ii), 0, 1);
    motor_noise(ii) = stdMotor;
end

figure
bar(motor_noise)
xlabel('Subject')
ylabel('Motor noise (deg)')