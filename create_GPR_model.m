function model = create_GPR_model(x_train,y_train,sigma,kfcn)
    % This function creates a structure. To predict on new points, use
    % y_new = model.pred(x_new)
    % To compute the gradient of the prediction at x_new with respect to
    % the training points, use 
    % d(y_new)/d(y_train) = model.grad(x_new)
    model.alpha = (kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)))\y_train;
    model.pred = @(x_pred) kfcn(x_pred,x_train)*model.alpha;
    model.grad = @(x_pred) kfcn(x_pred,x_train)*(kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)));
end