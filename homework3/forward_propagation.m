function activations = forward_propagation(X, parameters)
% Evelyn Kim, 706180341
% 
% forward_propagation: a function that computes the output 
% of a feedforward neural network given input data and learned parameters
% Inputs:
%       X: input data, shape (input size, number of examples)
%       parameters: learned parameters, a struct containing W1, b1, W2, b2, etc.
% Outputs:
%       activations: array of activations at each layer, including input and output layers

    L = length(parameters);
    A = X;
    activations = cell(1, L+1);
    activations{1} = X;
    for l = 1:L
        Z = parameters{l}.W*A + parameters{l}.b;
        if l == L
            A = softmax(Z);
        else
            A = tanh2(Z);
        end
        activations{l+1} = A;
    end
end