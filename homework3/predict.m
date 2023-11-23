function Y_pred = predict(X, parameters)
% Evelyn Kim, 706180341
% 
% predict: a function that returns the predicted classes of the inference images using 
% parameters containing the weights and biases of the neural network.
% Inputs:
%         X: Inference images with sizes 784 x N. N is the number of images. 
%         parameters: a struct containing the weights and biases (W1, b1,
%         W2, b2, etc.)
% Output: 
%         Y_pred: a 10 x N array containing the predicted labels of each input image.
    forward_pass = forward_propagation(X, parameters);
    Y_pred = forward_pass{end} == max(forward_pass{end});
end