function model = create_GPR_model2(x_train,y_train,sigma,kfcn,kfcn_grad)
    % This function creates a structure. To predict on new points, use
    % y_new = model.pred(x_new)
    % To compute the gradient of the prediction at x_new with respect to
    % the training points, use 
    % d(y_new)/d(y_train) = model.grad(x_new)
    model.kfcn = kfcn;
    model.alpha = (kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)))\y_train;
    model.pred = @(x_pred) kfcn(x_pred,x_train)*model.alpha;
    model.grad_dy = @(x_pred) kfcn(x_pred,x_train)/(kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)));
%     model.grad_dxdy = @(x_pred) kfcn_grad(x_pred,x_train)/(kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)));
    M = (kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)));
    model.grad_dxdy = @(x_pred) cat(3,kfcn_grad(x_pred,x_train,1)/M, kfcn_grad(x_pred,x_train,2)/M);
%     model.grad_dx = @(x_pred) kfcn_grad(x_pred,x_train)*model.alpha;
    model.grad_dx = @(x_pred) cat(2,kfcn_grad(x_pred,x_train,1)*model.alpha,kfcn_grad(x_pred,x_train,2)*model.alpha);
end