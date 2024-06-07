


#####function
generate_umap_data=function(mat,samp,clu){
       umap_genus=umap(t(log10(mat+0.000001)))
       umap_genus_data=data.frame(umap_genus$layout)
       umap_genus_data$ExpID3=rownames(umap_genus_data)
       umap_genus_data=merge(umap_genus_data,samp)
       dist_umap_genus=hclust(dist(umap_genus$layout),method="ward.D2")
       umap_genus_data$cluster=cutree(dist_umap_genus,clu)[umap_genus_data$ExpID3]
       return(umap_genus_data)
   }

hclusCut <- function(x, k, d.meth = "euclidean", h.meth="ward.D2",...)
  list(cluster = cutree(hclust(dist(x, method=d.meth), method=h.meth, ...), k=k))

fisher_matrix=function(mat){
  # p.adjust.method <- match.arg(p.adjust.method)
  dat=melt(mat)
  colnames(dat)[1:2]=c("V1","V2")
  dat$V1_sum=apply(mat,1,sum)[dat$V1]
  dat$V2_sum=apply(mat,2,sum)[dat$V2]
  dat$total=sum(dat$value)
  
  dat$V1_other=dat$V1_sum-dat$value
  dat$V2_other=dat$V2_sum-dat$value
  dat$V1_t_other=dat$total+dat$value-dat$V2_sum-dat$V1_sum
  
  dat$pvalue=apply(dat,1,function(a){fisher.test(matrix(as.numeric(a[c("value","V1_other","V2_other","V1_t_other")]),nrow = 2))$p.value})
  dat$estimate=apply(dat,1,function(a){fisher.test(matrix(as.numeric(a[c("value","V1_other","V2_other","V1_t_other")]),nrow = 2))$est})
  dat$qvalue=p.adjust(dat$pvalue,method = "fdr")
  return(dat)
}

setwd("E:/Project/2022_ECMO/ECMO_2024")
va_mat=read.csv("VA_case2.csv",row.names = "X")
va_case=read.csv("sample_va.csv",row.names = "X")
##sample all
dim(result_all_2023_genus_perc)


umap_genus_mat_va=generate_umap_data(va_mat,va_case,4)
ggplot(umap_genus_mat_va,aes(x=X1,y=X2,col=factor(cluster)))+geom_point()+theme_classic()

print_table(umap_genus_mat_va,"outcome","cluster")
fisher_matrix( table(umap_genus_mat_va[,c("outcome","cluster")]))
fisher.test(table(subset(umap_genus_mat_all_va,if_first=="Y")[,c("cluster.x","group_1219")]),simulate.p.value = T)

va_case=read.csv("va_case_p.csv")
seqsheet_blood_all_va=subset(seqsheet_blood_uniq_patient_res,patientID %in% va_case$patientID)
mat_genus_perc_all_va=result_all_2024_genus_perc[,seqsheet_blood_all_va$ExpID3]

umap_genus_mat_all_va=generate_umap_data(mat_genus_perc_all_va,seqsheet_blood_all_va,4)

dataset_va <- microtable$new(sample_table =umap_genus_mat_va,otu_table = va_mat, tax_table = otu_name_g_u)
lefse_va_umap <- trans_diff$new(dataset = dataset_va, 
                                method = "lefse", 
                                group = "cluster", 
                                alpha = 0.01, 
                                lefse_subgroup = NULL)

lefse_va_umap$plot_diff_cladogram(use_taxa_num = 200,use_feature_num = 50,clade_label_level = 6)

lefse_va_umap_genus <- trans_diff$new(dataset = dataset_va, 
  method = "lefse", 
  group = "cluster", 
  alpha = 0.01, 
  taxa_level = "genus")

seqsheet_blood_all_va=merge(va_case,seqsheet_blood_uniq_patient_res,by="ExpID3")
seqsheet_blood_all_va=subset(seqsheet_blood_uniq_patient_res,patientID %in% seqsheet_blood_all_va$patientID)
mat_genus_perc_all_va=result_all_2024_genus_perc[,seqsheet_blood_all_va$ExpID3]
umap_genus_mat_all_va=generate_umap_data(mat_genus_perc_all_va,seqsheet_blood_all_va,4)


first_day=data.frame(days=tapply(subset(umap_genus_mat_all_va,label =="E")$days,subset(umap_genus_mat_all_va,label =="E")$patientID,min))
first_day$patientID=rownames(first_day)
first_day$label="E"
first_day$if_first="Y"
umap_genus_mat_all_va=merge(umap_genus_mat_all_va,first_day,all.x = T)
umap_genus_mat_all_va=merge(umap_genus_mat_all_va,subset(umap_genus_mat_all_va,if_first=="Y")[,c("patientID","cluster")],by="patientID")
umap_genus_mat_all_va$cluster.x=factor(umap_genus_mat_all_va$cluster.x,levels = c(1,4,2,3),labels = c("A","B","C","D"))
umap_genus_mat_all_va$cluster.y=factor(umap_genus_mat_all_va$cluster.y,levels = c(1,4,2,3),labels = c("A","B","C","D"))

p1=ggplot(umap_genus_mat_all_va,aes(x=X1,y=X2,col=cluster.x))+geom_point()+theme_classic()
p2=ggplot(subset(umap_genus_mat_all_va,if_first=="Y"),aes(x=X1,y=X2,col=cluster.x))+geom_point()+theme_classic()
p1+p2

umap_genus_mat_all_va[umap_genus_mat_all_va$days==0,"days_new"]=1
umap_genus_mat_all_va[umap_genus_mat_all_va$days>7,"days_new"]=7

patient_days=data.frame(E_min=tapply(subset(umap_genus_mat_all_va,label =="E")$days_new,subset(umap_genus_mat_all_va,label =="E")$patientID,min))
patient_days$patientID=rownames(patient_days)
patient_days$E_max=tapply(subset(umap_genus_mat_all_va,label =="E")$days_new,subset(umap_genus_mat_all_va,label =="E")$patientID,max)[patient_days$patientID]
patient_days$label="E"
patient_days$P_min=tapply(subset(umap_genus_mat_all_va,label =="P")$days_new,subset(umap_genus_mat_all_va,label =="P")$patientID,min)[patient_days$patientID]
patient_days$P_max=tapply(subset(umap_genus_mat_all_va,label =="P")$days_new,subset(umap_genus_mat_all_va,label =="P")$patientID,max)[patient_days$patientID]
patient_days$label2="P"
patient_days=merge(patient_days,va_case[,c("patientID","survival","outcome")])
patient_days=merge(patient_days,subset(umap_genus_mat_all_va,if_first=="Y")[,c("patientID","cluster.x")])

patient_days=patient_days[order(patient_days$P_max ),]
patient_days=patient_days[order(patient_days$E_max ),]
patient_days=patient_days[order(patient_days$cluster.x),]

patient_days_before=patient_days[,c(1:4,8:10)]
patient_days_after=patient_days[,c(1,5:7,8:10)]
colnames(patient_days_before)[2:4]=c("min","max","label")
colnames(patient_days_after)[2:4]=c("min","max","label")

patient_days_l=rbind(patient_days_before,patient_days_after)

ggplot()+geom_segment(data=patient_days_l,aes(x=min,xend=max,y=patientID,yend=patientID),col="grey")+geom_point(data=umap_genus_mat_all_va,aes(x=days_new,y=patientID,col=factor(cluster.x)))+facet_grid(~label ,scales = "free",space = "free")+geom_point(data=subset(patient_days_l,(outcome=="die_before_wean"&label=="E")|(outcome!="die_before_wean"&label=="P")),aes(x=max,y=patientID,shape=outcome))+theme_classic()+ylim(rev(patient_days$patientID))+scale_color_brewer(palette = "Set1")

rownames(umap_genus_mat_all_va)=umap_genus_mat_all_va$ExpID3
dataset_va <- microtable$new(sample_table =subset(umap_genus_mat_all_va,if_first=="Y"),otu_table = mat_genus_perc_all_va[,subset(umap_genus_mat_all_va,if_first=="Y")$ExpID3], tax_table = otu_name_g_u)
lefse_va_umap <- trans_diff$new(dataset = dataset_va, 
                                 method = "lefse", 
                                 group = "cluster.x", 
                                 alpha = 0.01, 
                                 lefse_subgroup = NULL)

mt_genus_va_pa=melt(mat_genus_perc_all_va[,subset(umap_genus_mat_all_va,if_first=="Y")$ExpID3])
mt_genus_va_pa=subset(mt_genus_va_pa,value>0)
mt_genus_va_pa=merge(mt_genus_va_pa,subset(umap_genus_mat_all_va,if_first=="Y"),by.x="seqID",by.y="ExpID3")

genus_mean=tapply(mt_genus_va_pa$value,mt_genus_va_pa$short_genus,sum)
genus_mean=genus_mean[order(genus_mean,decreasing = T)]
pheatmap(mat_genus_perc_all_va[rownames(head(genus_mean,20)),subset(umap_genus_mat_all_va,if_first=="Y")$ExpID3],annotation_col =subset(umap_genus_mat_all_va,if_first=="Y")[,c("group_1219","cluster.x")],cluster_cols = F,show_colnames = F)

lefse_va_umap_genus_fir <- trans_diff$new(dataset = dataset_va, 
                                method = "lefse", 
                                group = "cluster.x", 
                                alpha = 0.01, taxa_level = "Genus",
                                lefse_subgroup = NULL)


calculate_beta=function(va_mat,va_case,group_select){
 library(vegan)
library(ape) 
beta_table <- as.matrix(vegdist(t(va_mat)), method = "bray", na.rm = F)
PCOA <- pcoa(beta_table)$vectors
var_exp <- pcoa(beta_table)$values
# Run stats for differentiation centroids
beta_dist = as.dist(beta_table)
length(beta_dist)
# Run PERMANOVA
ad = adonis(beta_dist ~ va_case[,group_select], permutations=999)

return(ad$aov.tab)
}


###WGS data
kraken_blood=read.csv("kraken_blood.csv",row.names = "X")
kraken_blood$perc=kraken_blood$readsCount/kraken_blood$all_kraken_bac_reads
genus_meta_count=tapply(kraken_blood$perc,kraken_blood$genus,sum)
genus_meta_count=genus_meta_count[order(genus_meta_count,decreasing = T)]
ggplot(subset(kraken_blood,genus  %in% rownames(head(genus_meta_count,10))),aes(y= EID ,x=perc,fill=genus))+geom_bar(stat="identity")+theme_classic()



mat_genus_exp=tapply(kraken_blood$readsCount,kraken_blood[,c("genus","EID")],sum)
mat_genus_exp[is.na(mat_genus_exp)]=0
mat_genus_exp_perc=t(t(mat_genus_exp)/apply(mat_genus_exp,2,sum))
mat_genus_exp_perc_top10=mat_genus_exp_perc[rownames(head(genus_meta_count,10)),]
mat_genus_exp_perc_top10=rbind(mat_genus_exp_perc_top10,1-apply(mat_genus_exp_perc_top10,2,sum))
rownames(mat_genus_exp_perc_top10)[11]="Other"

mat_genus_exp_perc_top10=t(mat_genus_exp_perc_top10)
mat_genus_exp_perc_top10=mat_genus_exp_perc_top10[order(mat_genus_exp_perc_top10[,6],decreasing = T),]
mat_genus_exp_perc_top10=mat_genus_exp_perc_top10[order(mat_genus_exp_perc_top10[,5],decreasing = T),]
mat_genus_exp_perc_top10=mat_genus_exp_perc_top10[order(mat_genus_exp_perc_top10[,4],decreasing = T),]
mat_genus_exp_perc_top10=mat_genus_exp_perc_top10[order(mat_genus_exp_perc_top10[,3],decreasing = T),]
mat_genus_exp_perc_top10=mat_genus_exp_perc_top10[order(mat_genus_exp_perc_top10[,2],decreasing = T),]
mat_genus_exp_perc_top10=mat_genus_exp_perc_top10[order(mat_genus_exp_perc_top10[,1],decreasing = T),]
ggplot(subset(melt(mat_genus_exp_perc_top10),value>0),aes(y= Var1 ,x=value,fill=Var2))+geom_bar(stat="identity")+theme_classic()+scale_fill_brewer(palette = "Set3")

metaphlan_blood=read.csv("metaphlan_blood.csv",row.names = "X")
ggplot(metaphlan_blood,aes(y= ExpID ,x=perc,fill=genus))+geom_bar(stat="identity")+theme_classic()


sample_wgs_16s=read.csv("E:/Project/202402_skw/blood_sample_16S_f.csv")
sample_wgs_16s=merge(sample_wgs_16s,umap_genus_mat,by.x="ExpID_16S",by.y="ExpID3",all.x = T)
rownames(sample_wgs_16s)=sample_wgs_16s$ExpID
kraken_wgs_mat=tapply(kraken_blood$perc, kraken_blood[,c("genus","EID")], sum)
pheatmap(kraken_wgs_mat[apply(kraken_wgs_mat,1,sum)>0.1,],annotation_col = sample_wgs_16s[,c("outcome","cluster")])


###add new batch

mat_16s_2402=read.delim("16S24DTH2_all.out",row.names = "Descrpiton")
mat_16s_2407=read.delim("16S24DTH7_all.out",row.names = "Descrpiton")

mt_16s_2402=subset(melt(as.matrix(mat_16s_2402)),value>0)
mt_16s_2407=subset(melt(as.matrix(mat_16s_2407)),value>0)

write.csv(rbind(mt_16s_2402,mt_16s_2407),file="mt_16s_long.csv")

result_long_16S24=read.csv("mt_24_long.csv")
result_all_2024=rbind(result_all_2023,result_long_16S24)
result_all_2024_g=subset(result_all_2024,Genus!="g__"&Genus!="")
otu_name_g_u=unique(result_all_2023_g[,c(4:9,11)])
rownames(otu_name_g_u)=otu_name_g_u$short_genus
seqsheet_blood_uniq_patient_res=read.csv("seqsheet_blood_uniq_patient_res.csv",row.names = "X")

result_all_2024_g_rm=subset(result_all_2024_g,!short_genus%in%c("f__Vibrionaceae;g__Salinivibrio"))

result_all_2024_genus=tapply(result_all_2024_g_rm$count,result_all_2024_g_rm[,c("short_genus","seqID")],sum)
result_all_2024_genus[is.na(result_all_2024_genus)]=0

result_all_2024_genus_perc=t(t(result_all_2024_genus)/apply(result_all_2024_genus,2,sum))

mat_genus_all=result_all_2024_genus[,seqsheet_blood_uniq_patient_res$ExpID3]
mat_genus_perc_all=result_all_2024_genus_perc[,seqsheet_blood_uniq_patient_res$ExpID3]

umap_genus_mat=generate_umap_data(mat_genus_perc_all,seqsheet_blood_uniq_patient_res,4)
p1=ggplot(umap_genus_mat,aes(x=X1,y=X2,col=factor(lib)))+geom_point()+theme_classic()+scale_color_brewer(palette = "Set1")+ggtitle("UMAP log")

beta_table_all <- as.matrix(vegdist(t(mat_genus_perc_all)), method = "bray", na.rm = F)
PCOA_all <- pcoa(beta_table_all)$vectors
pcoa_data_all=data.frame(PCOA_all[,1:2])
pcoa_data_all$ExpID3=rownames(pcoa_data_all)
pcoa_data_all=merge(pcoa_data_all,seqsheet_blood_uniq_patient_res)
p2=ggplot(pcoa_data_all,aes(x=Axis.1,y= Axis.2,col=lib))+geom_point()+theme_classic()+scale_color_brewer(palette = "Set1")+ggtitle("PCOA Bary")

pca_genus=prcomp(t(mat_genus_perc_all))
pca_genus_data=data.frame(pca_genus$x[,1:2])
pca_genus_data$ExpID3=rownames(pca_genus_data)
pca_genus_data=merge(pca_genus_data,seqsheet_blood_uniq_patient_res)
p3=ggplot(pca_genus_data,aes(x=PC1,y= PC2,col=lib))+geom_point()+theme_classic()+scale_color_brewer(palette = "Set1")+ggtitle("PCA")

p1+p2+p3

sum_genus_count=apply(mat_genus_perc_all, 1, sum)
sum_genus_count=sum_genus_count[order(sum_genus_count,decreasing = T)]
pheatmap(mat_genus_perc_all[names(head(sum_genus_count,30)),],annotation_col = seqsheet_blood_uniq_patient_res[,c("group_1219","lib")],clustering_method = "ward.D2",show_colnames = F)

sum_genus_count_log=apply(log10(mat_genus_perc_all+0.000001), 1, sum)
sum_genus_count_log=sum_genus_count_log[order(sum_genus_count_log,decreasing = T)]
pheatmap(log10(mat_genus_perc_all+0.000001)[names(head(sum_genus_count_log,30)),],annotation_col = seqsheet_blood_uniq_patient_res[,c("group_1219","lib")],clustering_method = "ward.D2",show_colnames = F)


batch1_sample=subset(seqsheet_blood_uniq_patient_res,lib %in% c("ecmo_1","ecmo_2","ecmo_3","ecmo_4"))
batch3_sample=subset(seqsheet_blood_uniq_patient_res,lib %in% c("ecmo_12"))
batch1_mat=mat_genus_perc_all[,batch1_sample$ExpID3]
batch3_mat=mat_genus_perc_all[,batch3_sample$ExpID3]

tbl_beta_batch1=calculate_beta(batch1_mat,batch1_sample,"group_1219")
tbl_beta_batch3=calculate_beta(batch3_mat,batch3_sample,"group_1219")

PCOA_batch1 <- pcoa(as.matrix(vegdist(t(batch1_mat)), method = "bray", na.rm = F))$vectors
pcoa_data_batch1=data.frame(PCOA_batch1[,1:2])
pcoa_data_batch1$ExpID3=rownames(pcoa_data_batch1)
pcoa_data_batch1=merge(pcoa_data_batch1,batch1_sample)
p1=ggplot(pcoa_data_batch1,aes(x=Axis.1,y= Axis.2,col=group_1219 ))+geom_point()+stat_ellipse(level = 0.68)+theme_classic()+scale_color_brewer(palette = "Set1")

PCOA_batch3 <- pcoa(as.matrix(vegdist(t(batch3_mat)), method = "bray", na.rm = F))$vectors
pcoa_data_batch3=data.frame(PCOA_batch3[,1:2])
pcoa_data_batch3$ExpID3=rownames(pcoa_data_batch3)
pcoa_data_batch3=merge(pcoa_data_batch3,batch3_sample)
p2=ggplot(pcoa_data_batch3,aes(x=Axis.1,y= Axis.2,col=group_1219 ))+geom_point()+stat_ellipse(level = 0.68)+theme_classic()+scale_color_brewer(palette = "Set1")

batch1_sample_rm=subset(seqsheet_blood_uniq_patient_res,lib %in% c("ecmo_1","ecmo_2","ecmo_3","ecmo_4")&group_1219 !="die_before_wean")
batch3_sample_rm=subset(seqsheet_blood_uniq_patient_res,lib %in% c("ecmo_12")&group_1219 !="die_before_wean")
batch1_mat_rm=mat_genus_perc_all[,batch1_sample_rm$ExpID3]
batch3_mat_rm=mat_genus_perc_all[,batch3_sample_rm$ExpID3]

dataset_batch1_rm <- microtable$new(sample_table =batch1_sample_rm,otu_table = batch1_mat_rm, tax_table = otu_name_g_u)
lefse_batch1_rm <- trans_diff$new(dataset = dataset_batch1_rm, 
                                method = "lefse", 
                                group = "group_1219", 
                                alpha = 0.05, p_adjust_method = NULL,
                                lefse_subgroup = NULL,taxa_level = "Genus")

dataset_batch3_rm <- microtable$new(sample_table =batch3_sample_rm,otu_table = batch3_mat_rm, tax_table = otu_name_g_u)
lefse_batch3_rm <- trans_diff$new(dataset = dataset_batch3_rm, 
                               method = "lefse", 
                               group = "group_1219", 
                               alpha = 0.05, p_adjust_method = NULL,
                               lefse_subgroup = NULL,taxa_level = "Genus")

lefse_batch1_rm$res_diff[intersect(lefse_batch1_rm$res_diff$Taxa,lefse_batch3_rm$res_diff$Taxa),]
lefse_batch3_rm$res_diff[intersect(lefse_batch1_rm$res_diff$Taxa,lefse_batch3_rm$res_diff$Taxa),]

######cluster and other result

fisher_matrix( table(subset(umap_genus_mat_all_va,if_first=="Y")[,c("group_1219","cluster.x")]))
cyto_result=read.csv("cyto_result.csv")
cyto_result=merge(cyto_result,subset(umap_genus_mat_all_va,if_first=="Y")[,c("patientID_n","group_1219","cluster.x")],by.x="patientID",by.y="patientID_n",all.x = T)
colnames(cyto_result)[36]="cluster_16S"
cyto_table1=CreateTableOne(vars =  colnames(cyto_result)[c(6:9,25:34)],strata = "cluster_16S",data=subset(cyto_result,!is.na(cluster_16S)))
tbl_cyto=print(cyto_table1,nonnormal = T )
write.csv(tbl_cyto,"tbl_cyto.csv")

ecmo_cyto_all=read.csv("ecmo_cyto.csv",row.names = "X")
ecmo_cyto_all=merge(ecmo_cyto_all[,c(1:25)],subset(umap_genus_mat_all_va,if_first=="Y")[,c("patientID_n","group_1219","cluster.x")],by.x="patientID",by.y="patientID_n",all.x = T)
colnames(ecmo_cyto_all)[27]="cluster_16S"
cyto_all_table1=CreateTableOne(vars = colnames(ecmo_cyto_all)[7:25],strata = "cluster_16S",data=subset(ecmo_cyto_all,!is.na(cluster_16S)))
tbl_cyto_all=print(cyto_all_table1,nonnormal = T )
write.csv(tbl_cyto_all,"tbl_cyto_all.csv")

cd4_cluster_result=merge(cd4_cluster_result,subset(umap_genus_mat_all_va,if_first=="Y")[,c("patientID_n","group_1219","cluster.x")],by.x="patientID",by.y="patientID_n",all.x = T)
cd4_cluster_result=unique(cd4_cluster_result)
colnames(cd4_cluster_result)[29]="cluster_16S"
cd4_table1=CreateTableOne(vars = colnames(cd4_cluster_result)[2:27],strata = "cluster_16S",data=subset(cd4_cluster_result,!is.na(cluster_16S)))
tbl_cd4=print(cd4_table1,nonnormal = T )
write.csv(tbl_cd4,"tbl_cd4.csv")

elisa_result=read.csv("elisa_result.csv")
elisa_result=merge(elisa_result,subset(umap_genus_mat_all_va,if_first=="Y")[,c("patientID_n","group_1219","cluster.x")],by.x="sample",by.y="patientID_n",all.x = T)
colnames(elisa_result)[9]="cluster_16S"
elisa_table1=CreateTableOne(vars = colnames(elisa_result)[4:7],strata = "cluster_16S",data=subset(elisa_result,!is.na(cluster_16S)))
tbl_elisa=print(elisa_table1,nonnormal = T )
write.csv(tbl_elisa,"tbl_elisa.csv")

ggplot(subset(elisa_result,!is.na(cluster_16S)),aes(x=cluster_16S,y=LBP,fill= cluster_16S))+geom_boxplot()+theme_classic()+stat_compare_means(comparisons = list(c("A","B"),c("A","C"),c("A","D"),c("B","C"),c("B","D"),c("C","D")),aes(label = ..p.signif..))
ggplot(subset(elisa_result,!is.na(cluster_16S)),aes(x=cluster_16S,y=FABP    ,fill= cluster_16S))+geom_boxplot()+theme_classic()+stat_compare_means(comparisons = list(c("A","B"),c("A","C"),c("A","D"),c("B","C"),c("B","D"),c("C","D")),aes(label = ..p.signif..))
ggplot(subset(elisa_result,!is.na(cluster_16S)),aes(x=cluster_16S,y=sCD14     ,fill= cluster_16S))+geom_boxplot()+theme_classic()+stat_compare_means(comparisons = list(c("A","B"),c("A","C"),c("A","D"),c("B","C"),c("B","D"),c("C","D")),aes(label = ..p.signif..))
ggplot(subset(elisa_result,!is.na(cluster_16S)),aes(x=cluster_16S,y=zonulin      ,fill= cluster_16S))+geom_boxplot()+theme_classic()+stat_compare_means(comparisons = list(c("A","B"),c("A","C"),c("A","D"),c("B","C"),c("B","D"),c("C","D")),aes(label = ..p.signif..))

