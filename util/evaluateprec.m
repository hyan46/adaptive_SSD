function [ fpr,fnr,precision,recall,tp,fp,tn,fn ] = evaluateprec(predIntry,trueIntry,nS)
%UNTITLED2 Summary of this function goes here
%%   Detailed explanation goes here

totalIntry = 1:nS;

% predIntryFalse1 = find(S == 0);
% trueIntryFalse1 = find(Strue == 0);

predIntryFalse = setdiff(totalIntry,predIntry);
trueIntryFalse = setdiff(totalIntry,trueIntry);

tp = length(intersect(predIntry,trueIntry));
fp = length(intersect(predIntry,trueIntryFalse));
fn = length(intersect(predIntryFalse,trueIntry));
tn = length(intersect(predIntryFalse,trueIntryFalse));


fpr = fp/(tn + fp);
fnr = fn/(tp + fn);
precision = tp/(tp + fp);
recall = tp/(tp + fn);
% accuracy = (tp + tn)/(tp + tn + fp + fn);





end

