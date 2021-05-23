function alpha = get_alpha(x_train,y_train,kfcn,sigma,train_format)
    %     if strcmp(train_format,'scattered')
    %         alpha = (kfcn(x_train,x_train,train_format) + sigma^2*eye(size(x_train,1),size(x_train,1)))\y_train;
    %     elseif strcmp(train_format,'gridded')
    %         alpha = (kfcn(x_train,x_train,train_format) + sigma^2*eye(size(x_train,1)^2,size(x_train,1)^2))\y_train;
    %     end
    alpha = (kfcn(x_train,x_train,train_format) + sigma^2*eye(size(x_train,1),size(x_train,1)))\y_train;
end