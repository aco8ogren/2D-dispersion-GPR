function grad_dxdy = get_grad_dxdy(x_pred,x_train,sigma,kfcn,kfcn_grad,query_format)
    %     if strcmp(query_format,'scattered')
    %         M = (kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)));
    %     elseif strcmp(query_format,'gridded')
    %         M = (kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1)^2,size(x_train,1)^2));
    %     end
    M = (kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)));
    grad_dxdy = cat(3,...
        kfcn_grad(x_pred,x_train,1)/M,...
        kfcn_grad(x_pred,x_train,2)/M);
end