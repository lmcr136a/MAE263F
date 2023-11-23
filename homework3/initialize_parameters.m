function parameters = initialize_parameters(layer_dims)
% Evelyn Kim, 706180341
% 
% initialize_parameters: a function that initializes the weights 
% and biases of the feedforward neural network
% Inputs:
%   layer_dims: array of layer dimensions, including input and output layers
% Output: 
%   parameters: a struct containing W1, b1, W2, b2, etc.

    L = length(layer_dims);
    parameters = cell(1, L-1);
    for l = 1:(L-1)
        parameters{l}.W  = randn(layer_dims(l+1), layer_dims(l));
        parameters{l}.b = zeros(layer_dims(l+1), 1);
    end
end