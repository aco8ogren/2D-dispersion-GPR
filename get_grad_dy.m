function grad_dy = get_grad_dy(x_pred,x_train,sigma,kfcn,query_format)
    %     if strcmp(query_format,'scattered')
    %         grad_dy = kfcn(x_pred,x_train,query_format)/(kfcn(x_train,x_train,query_format) + sigma^2*eye(size(x_train,1),size(x_train,1)));
    %     elseif strcmp(query_format,'gridded')
    %         grad_dy = kfcn(x_pred,x_train,query_format)/(kfcn(x_train,x_train,query_format) + sigma^2*eye(size(x_train,1)^2,size(x_train,1)^2));
    %     end
    grad_dy = kfcn(x_pred,x_train,query_format)/...
        (kfcn(x_train,x_train,query_format) + sigma^2*eye(size(x_train,1),size(x_train,1)));
end