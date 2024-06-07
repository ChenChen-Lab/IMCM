
data_sur=read.csv("data_sur.csv",row.names = "X")
lac_data=read.csv("completedData.csv")
lac_data$lac_z=(lac_data$lac-mean(lac_data$lac))/sd(lac_data$lac)

bootstrap_rf=function(data_sur,tbl_select,num_boot){
rf_model_predict=data.frame()

for(i in 1:num_boot){
  set.seed(i)
  sampled_rows <- sample(1:nrow(data_sur), size = 359*0.7, replace = FALSE)
  train_table=data_sur[sampled_rows, ]

  
  validataion_table=data_sur[-sampled_rows, ]
  
  
  features <- train_table[, tbl_select]  
  response <- factor(train_table$survival)  
  # 训练随机森林模型
  rf_model <- randomForest(x = features, y = response, importance = TRUE)
 
 tbl_train=table((data.frame(response,predict=predict(rf_model, features))))
 sens_train=tbl_train[1,1]/(tbl_train[1,1]+tbl_train[1,2])
 spec_train=tbl_train[2,2]/(tbl_train[2,1]+tbl_train[2,2])
 accu_train=(tbl_train[2,2]+tbl_train[1,1])/sum(tbl_train)
 
  features_validation=validataion_table[,tbl_select]
  response_validation=validataion_table$survival  
  rf_prob=predict(rf_model, features_validation,type="prob")
  roc_test=roc(response_validation,rf_prob[,2])
  auc_test=as.numeric(auc(roc_test))
  tbl_validation= table((data.frame(response_validation,predict=predict(rf_model, features_validation))))
  sens_val=tbl_validation[1,1]/(tbl_validation[1,1]+tbl_validation[1,2])
  spec_val=tbl_validation[2,2]/(tbl_validation[2,1]+tbl_validation[2,2])
  accu_val=(tbl_validation[2,2]+tbl_validation[1,1])/sum(tbl_validation)
  rf_model_predict=rbind(rf_model_predict,c(i,sens_train,spec_train,accu_train,sens_val,spec_val,accu_val,auc_test))
}
colnames(rf_model_predict)=c("seed","sens_train","spec_train","accu_train","sens_val","spec_val","accu_val","auc_val")
return(rf_model_predict)
}

rf_model_predict_bact=bootstrap_rf(data_sur,2:18,10)
data_sur_lac=merge(data_sur,lac_data[,c("ExpID3","lac")])
rf_model_predict_lac=bootstrap_rf(data_sur_lac,c(2:18,20),10)

lac_data$lac_z=(lac_data$lac-mean(lac_data$lac))/sd(lac_data$lac)
data_sur_lac_z=merge(data_sur,lac_data[,c("ExpID3","lac_z")])
rf_model_predict_lac_z=bootstrap_rf(data_sur_lac_z,c(2:18,20),100)

ggplot(rf_model_predict_lac_z)+geom_boxplot(aes(x="accuracy_train",y=accu_train  ))+geom_boxplot(aes(x="accuracy_validation",y=accu_val  ))+geom_jitter(aes(x="accuracy_train",y=accu_train  ))+geom_jitter(aes(x="accuracy_validation",y=accu_val  ))+geom_boxplot(aes(x="AUC_validation",y=auc_val ))+geom_jitter(aes(x="AUC_validation",y=auc_val  ))+theme_classic()
ggplot(rf_model_predict_lac_z)+geom_density(aes(x=auc_val ))+theme_classic()

set.seed(29)
data_sur=data_sur_lac_z
tbl_select=c(2:18,20)

sampled_rows <- sample(1:nrow(data_sur), size = 359*0.7, replace = FALSE)
train_table=data_sur[sampled_rows, ]


validataion_table=data_sur[-sampled_rows, ]


features <- train_table[, tbl_select]  
response <- factor(train_table$survival)  
# 训练随机森林模型
rf_model <- randomForest(x = features, y = response, importance = TRUE)

tbl_train=table((data.frame(response,predict=predict(rf_model, features))))
sens_train=tbl_train[1,1]/(tbl_train[1,1]+tbl_train[1,2])
spec_train=tbl_train[2,2]/(tbl_train[2,1]+tbl_train[2,2])
accu_train=(tbl_train[2,2]+tbl_train[1,1])/sum(tbl_train)

features_validation=validataion_table[,tbl_select]
response_validation=validataion_table$survival  
rf_prob=predict(rf_model, features_validation,type="prob")
roc_test=roc(response_validation,rf_prob[,2])
auc_test=as.numeric(auc(roc_test))

plot(roc_test)