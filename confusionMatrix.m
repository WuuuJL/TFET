function [TP, FN, FP, TN] = confusionMatrix(Label1, Label2, Prediction)
%   Input
%       Label1     --- label of impulsive-like components  -1
%       Label2     --- label of harmonic-like components    1
%       Prediction --- result of segmentation

%   Output
%       TP --- true positive
%       FN --- false negative
%       FP --- false positive
%       TN --- true negative
total_pixels = sum(sum(abs(Label1))) + sum(sum(Label2));

Tem1 = Prediction-Label1;
Tem2 = Prediction-Label2;
FN = sum(sum(Tem1 == (2+zeros(size(Tem1)))));
FP = sum(sum(Tem2 == (-2+zeros(size(Tem2)))));
TP = sum(sum((Prediction == Label1) & (Label1 == -1)));
TN = sum(sum((Prediction == Label2) & (Label2 == 1)));

if TP+TN+FN+FP ~= total_pixels
    disp('There is an error ');
end
end