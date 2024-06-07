kraken_ases=read.csv("kraken_sp_ascites.csv")
select_genus=c("g__Klebsiella","g__Pseudomonas","g__Escherichia","g__Enterococcus","g__Staphylococcus")
kraken_ases_genus=subset(kraken_ases,genus %in% select_genus)

kraken_genus_count_mt=melt(tapply(kraken_ases_genus$readsCount,kraken_ases_genus[,c("ExpID","genus")],sum))
kraken_genus_count_mt[is.na(kraken_genus_count_mt$value),"value"]=0
#kraken_genus_count_mt=subset(kraken_genus_count_mt,!is.na(value))
kraken_genus_count_mt$sum=tapply(kraken_ases_genus$readsCount,kraken_ases_genus[,c("ExpID")],sum)[kraken_genus_count_mt$ExpID]
kraken_genus_count_mt$perc=kraken_genus_count_mt$value/kraken_genus_count_mt$sum

ases_culture=read.csv("ases_culture.csv")


kraken_genus_count_mt_culture_pos=subset(kraken_genus_count_mt,ExpID %in% ases_culture$ExpID)
kraken_genus_count_mt_culture_pos=merge(kraken_genus_count_mt_culture_pos,ases_culture[,c(1,3,4)],all.x = T)
kraken_genus_count_mt_culture_pos[is.na(kraken_genus_count_mt_culture_pos$culture),"culture"]="N"



data_sene_genus=function(kraken_genus_count_mt_culture_pos,select_genus){
  library(pROC)
coords_sg=data.frame()
data_sg=data.frame()
for (sg in select_genus){
  roc_sg=roc(subset(kraken_genus_count_mt_culture_pos,genus==sg)$culture,subset(kraken_genus_count_mt_culture_pos,genus==sg)$perc)
  add_data_sg=data.frame(coords(roc_sg))
  add_data_sg$genus=sg
  data_sg=rbind(data_sg,add_data_sg)
  coords_sg=rbind(coords_sg,c(as.numeric(auc(roc_sg)), as.numeric(coords(roc_sg, "best", ret=c("threshold", "specificity", "sensitivity"))[1,]),sg))
}
colnames(coords_sg)=c("auc","threshold","specificity", "sensitivity","genus")
coords_sg$sensitivity=as.numeric(coords_sg$sensitivity)
coords_sg$specificity=as.numeric(coords_sg$specificity)
return(coords_sg)
}

ggplot()+geom_line(data=data_sg,aes(y=sensitivity ,x=specificity ,col=genus))+theme_classic()+geom_point(data=coords_sg,aes(y=sensitivity ,x=specificity ,col=genus))+geom_text(data=coords_sg,aes(y=sensitivity ,x=specificity ,label=paste(genus,threshold ,sep="\n")))


kraken_genus_count_mt_culture_pos=merge(kraken_genus_count_mt_culture_pos,coords_sg[,c("genus","threshold")])
kraken_genus_count_mt_culture_pos$ifPos="N"
kraken_genus_count_mt_culture_pos[kraken_genus_count_mt_culture_pos$perc>kraken_genus_count_mt_culture_pos$threshold,"ifPos"]="Y"

kraken_genus_count_mt=merge(kraken_genus_count_mt,coords_sg[,c("genus","threshold")])
kraken_genus_count_mt$ifPos="N"
kraken_genus_count_mt[kraken_genus_count_mt$perc>kraken_genus_count_mt$threshold,"ifPos"]="Y"


blood_culture=read.csv("blood_culture_result.csv")
colnames(blood_culture)[1:2]=c("Patient","genus")
blood_kraken=read.csv("kraken_blood.csv",row.names = "X")
colnames(blood_kraken)[1:2]=c("ExpID","Patient")

kraken_genus_count_mt=melt(tapply(blood_kraken$readsCount,blood_kraken[,c("Patient","genus")],sum))
kraken_genus_count_mt[is.na(kraken_genus_count_mt$value),"value"]=0
#kraken_genus_count_mt=subset(kraken_genus_count_mt,!is.na(value))
kraken_genus_count_mt$sum=tapply(blood_kraken$readsCount,blood_kraken[,c("Patient")],sum)[kraken_genus_count_mt$Patient]
kraken_genus_count_mt$perc=kraken_genus_count_mt$value/kraken_genus_count_mt$sum

blood_kraken_selected=subset(kraken_genus_count_mt,Patient%in% blood_culture$Patient & genus %in% blood_culture$genus)
blood_kraken_selected=merge(blood_kraken_selected,blood_culture[,c(1,2,4)],all.x="T")
blood_kraken_selected[is.na(blood_kraken_selected$culture),"culture"]="N"
blood_kraken_selected$genus=as.character(blood_kraken_selected$genus)
tbl_blood=table(blood_kraken_selected[,c("genus","culture")])
select_genus_bl=rownames(tbl_blood[tbl_blood[,2]>0,])
blood_kraken_selected=subset(blood_kraken_selected,genus %in% select_genus_bl)

coords_sg=data.frame()
data_sg=data.frame()
for (sg in select_genus_bl){
  roc_sg=roc(subset(blood_kraken_selected,genus==sg)$culture,subset(blood_kraken_selected,genus==sg)$perc)
  add_data_sg=data.frame(coords(roc_sg))
  add_data_sg$genus=sg
  data_sg=rbind(data_sg,add_data_sg)
  coords_sg=rbind(coords_sg,c(as.numeric(auc(roc_sg)), as.numeric(coords(roc_sg, "best", ret=c("threshold", "specificity", "sensitivity"))[1,]),sg))
}
colnames(coords_sg)=c("auc","threshold","specificity", "sensitivity","genus")
coords_sg$sensitivity=as.numeric(coords_sg$sensitivity)
coords_sg$specificity=as.numeric(coords_sg$specificity)

ggplot()+geom_line(data=data_sg,aes(y=sensitivity ,x=specificity ,col=genus))+theme_classic()+geom_point(data=coords_sg,aes(y=sensitivity ,x=specificity ,col=genus))+geom_text(data=coords_sg,aes(y=sensitivity ,x=specificity ,label=paste(genus,threshold ,sep="\n")))

blood_kraken_selected=merge(blood_kraken_selected,subset(coords_sg,auc>0.5)[,c("genus","threshold")])
blood_kraken_selected$ngs="N"
blood_kraken_selected[blood_kraken_selected$perc>blood_kraken_selected$threshold,"ngs"]="Y"