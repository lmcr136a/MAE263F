function cost = compute_cost(AL, Y)
% Evelyn Kim, 706180341
% 
% computer_cose: a function that calculates a cross-entropy loss between
% the model's prediction and the ground truth.
% Inputs:
%       AL: model's output. data shape: N*C, N is batch size and C is class number
%       Y: ground truth labels. data shape: N*C, N is batch size and C is class number
% Outputs:
%       returns: a constant value that sum of the loss of all inputs.
    cost = -sum(Y .* log(AL));
end