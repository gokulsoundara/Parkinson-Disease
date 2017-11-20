##Read the data in data frames for chromosome 1
cc_map <-read.table("C:/Users/gokul/Documents/Projects/Parkinson Disease/cc_map/chr11.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("C:/Users/gokul/Documents/Projects/Parkinson Disease/cc_pre/cc_chr11.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("C:/Users/gokul/Documents/Projects/Parkinson Disease/pd_map/chr11.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("C:/Users/gokul/Documents/Projects/Parkinson Disease/pd_pre/chr11.pre", sep=" ",stringsAsFactors = F)
#C:\Users\gokul\Documents\Projects\Parkinson Disease

##initize the counters
A=0
a=0
#aa=0
#ds
#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      A_cc a_cc  A_pd a_pd  chi_sq
#1           2 rs4637157   232   36     3   223    NA     
#3           2  rs300789   173   92     6   180    NA     
#6           2  rs300773   190   74     7   178    NA     
#8           2  rs300711   181   84     6   200    NA     
#9           2  rs381726    81  130    60    91    NA    
#10          2  rs408209    83  128    60    91    NA    

df <- data.frame(matrix(ncol = 7, nrow = 0))
x <- c("chromosome", "locus", "A_cc","a_cc","A_pd","a_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  A = sum(cc_pre[,2*i+2] == cc_map[i,4])
  a = sum(cc_pre[,2*i+2] == cc_map[i,5])
  A = A + sum(cc_pre[,2*i+3] == cc_map[i,4])
  a = a + sum(cc_pre[,2*i+3] == cc_map[i,5])
  df[i,"A_cc"] = A
  df[i,"a_cc"] = a
}


A=0
a=0
#Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  A = sum(pd_pre[,2*i+2] == pd_map[i,4])
  a = sum(pd_pre[,2*i+2] == pd_map[i,5])
  A = A + sum(pd_pre[,2*i+3] == pd_map[i,4])
  a = a + sum(pd_pre[,2*i+3] == pd_map[i,5])
  df[i,"A_pd"] = A
  df[i,"a_pd"] = a
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chr..omosome 1
df_combined=df
##################


###Row bind the dataframes for a cpmplete dataset for all chromosomes for all locii
#df_combined= rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22,dfX,dfY,dfXY)

############################
############################
##first define the function to apply
Chsq <- function(x){
  ## input is a row of your data
  ## creating a table from each row
  if (x[1]!=0 & x[2]=0 )  y <- matrix(x[c(1,3)],byrow =TRUE,nrow=2)
  if (x[1]=0 & x[2]!=0 )  y <- matrix(x[c(2,4)],byrow =TRUE,nrow=2)
  if (x[1]==0 & x[2]==0)  return (1) ##for all control zero count 
  
  ### this will return the p value
  return(chisq.test(y)$p.value)
}
data =df_combined[,3:6]
P_Values <- as.vector(apply(data,1,Chsq))   ###return p values
#df_combined['chi_sq'] = P_Values

df_combined = cbind(df_combined,P_Values) ####Populate p values in the df 

head(df_combined)
df_combined[df_combined$locus == "rs3741411", ]
y=matrix(c(442,100,261,9),byrow = TRUE,nrow = 2)
y