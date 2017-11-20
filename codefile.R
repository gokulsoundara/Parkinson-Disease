##Read the data in data frames for chromosome 1
cc_map <-read.table("C:/Users/gokul/Documents/Parkinson Disease/cc_map/chr1.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("C:/Users/gokul/Documents/Parkinson Disease/cc_pre/cc_chr1.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("C:/Users/gokul/Documents/Parkinson Disease/pd_map/chr1.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("C:/Users/gokul/Documents/Parkinson Disease/pd_pre/chr1.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0
ds
#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chr..omosome 1
df1=df
################################################################################
##Read the data in data frames for chromosome 2

cc_map <-read.table("./cc_map/chr2.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr2.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr2.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr2.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)
  
  for (i in 1:nrow(cc_map)) {
      
      Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
      AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
      aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
      df[i,"AA_cc"]=AA
      df[i,"Aa_cc"]=Aa
      df[i,"aa_cc"]=aa
    }
  

AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
      Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
      AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
      aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
    
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
 #populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]


##Store the dataframe df for chromosome 2
df2=df

###########################################################

##Read the data in data frames for chromosome 3
cc_map <-read.table("./cc_map/chr3.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr3.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr3.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr3.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 3
df3=df
############################################################

###########################################################

##Read the data in data frames for chromosome 4
cc_map <-read.table("./cc_map/chr4.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr4.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr4.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr4.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 3
df4=df
############################################################

###########################################################

##Read the data in data frames for chromosome 5
cc_map <-read.table("./cc_map/chr5.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr5.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr5.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr5.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 5
df5=df
############################################################

###########################################################

##Read the data in data frames for chromosome 1
cc_map <-read.table("./cc_map/chr6.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr6.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr6.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr6.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 6
df6=df
############################################################
###########################################################

##Read the data in data frames for chromosome 7
cc_map <-read.table("./cc_map/chr7.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr7.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr7.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr7.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 7
df7=df
############################################################
###########################################################

##Read the data in data frames for chromosome 8
cc_map <-read.table("./cc_map/chr8.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr8.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr8.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr8.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 8
df8=df
############################################################
###########################################################

##Read the data in data frames for chromosome 9
cc_map <-read.table("./cc_map/chr9.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr9.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr9.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr9.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 9
df9=df
############################################################
###########################################################

##Read the data in data frames for chromosome 10
cc_map <-read.table("./cc_map/chr10.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr10.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr10.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr10.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 10
df10=df
############################################################
###########################################################

##Read the data in data frames for chromosome 11
cc_map <-read.table("./cc_map/chr11.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr11.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr11.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr11.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 11
df11=df
############################################################
###########################################################

##Read the data in data frames for chromosome 12
cc_map <-read.table("./cc_map/chr12.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr12.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr12.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr12.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 12
df12=df
############################################################
###########################################################

##Read the data in data frames for chromosome 13
cc_map <-read.table("./cc_map/chr13.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr13.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr13.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr13.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 13
df13=df
############################################################
###########################################################

##Read the data in data frames for chromosome 14
cc_map <-read.table("./cc_map/chr14.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr14.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr14.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr14.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]


##Store the dataframe df for chromosome 14
df14=df
############################################################
###########################################################

##Read the data in data frames for chromosome 15
cc_map <-read.table("./cc_map/chr15.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr15.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr15.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr15.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 15
df15=df
############################################################
###########################################################

##Read the data in data frames for chromosome 16
cc_map <-read.table("./cc_map/chr16.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr16.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr16.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr16.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 16
df16=df
############################################################
###########################################################

##Read the data in data frames for chromosome 17
cc_map <-read.table("./cc_map/chr17.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr17.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr17.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr17.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 17
df17=df
############################################################
###########################################################

##Read the data in data frames for chromosome 18
cc_map <-read.table("./cc_map/chr18.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr18.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr18.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr18.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 18
df18=df
############################################################
###########################################################

##Read the data in data frames for chromosome 19
cc_map <-read.table("./cc_map/chr19.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr19.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr19.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr19.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 19
df19=df
############################################################
###########################################################

##Read the data in data frames for chromosome 20
cc_map <-read.table("./cc_map/chr20.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr20.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr20.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr20.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 20
df20=df
############################################################
###########################################################

##Read the data in data frames for chromosome 21
cc_map <-read.table("./cc_map/chr21.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr21.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr21.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr21.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]


##Store the dataframe df for chromosome 21
df21=df
############################################################
###########################################################

##Read the data in data frames for chromosome 22
cc_map <-read.table("./cc_map/chr22.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chr22.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chr22.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chr22.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome 22
df22=df
############################################################
###########################################################

##Read the data in data frames for chromosome X
cc_map <-read.table("./cc_map/chrX.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chrX.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chrX.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chrX.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]


##Store the dataframe df for chromosome X
dfX=df
############################################################
###########################################################

##Read the data in data frames for chromosome Y
cc_map <-read.table("./cc_map/chrY.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chrY.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chrY.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chrY.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]



##Store the dataframe df for chromosome Y
dfY=df
############################################################

###########################################################

##Read the data in data frames for chromosome XY
cc_map <-read.table("./cc_map/chrXY.map", sep="\t",stringsAsFactors = F)
cc_pre <-read.table("./cc_pre/cc_chrXY.pre", sep=" ",stringsAsFactors = F)
pd_map <-read.table("./pd_map/chrXY.map", sep="\t",stringsAsFactors = F)
pd_pre <-read.table("./pd_pre/chrXY.pre", sep=" ",stringsAsFactors = F)

##initize the counters
Aa=0
AA=0
aa=0

#create a data frame that will hold data based on chromosomes and locus in the follwing format:
#     chromosomelocus      AA_cc Aa_cc aa_cc AA_pd Aa_pd aa_pd chi_sq
#1           2 rs4637157   232    36     3   223    44     3     NA
#3           2  rs300789   173    92     6   180    81     9     NA
#6           2  rs300773   190    74     7   178    84     8     NA
#8           2  rs300711   181    84     6   200    63     7     NA
#9           2  rs381726    81   130    60    91   128    51     NA
#10          2  rs408209    83   128    60    91   131    48     NA

df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("chromosome", "locus", "AA_cc","Aa_cc","aa_cc","AA_pd","Aa_pd","aa_pd","chi_sq")
colnames(df) <- x

## populate in df the count of AA,Aa, and aa for control(cc)

for (i in 1:nrow(cc_map)) {
  
  Aa=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  AA=sum(cc_map[i,4]==cc_pre[,2*i+2] & cc_map[i,4]==cc_pre[,2*i+3]) 
  aa=sum(cc_map[i,5]==cc_pre[,2*i+2] & cc_map[i,5]==cc_pre[,2*i+3]) 
  df[i,"AA_cc"]=AA
  df[i,"Aa_cc"]=Aa
  df[i,"aa_cc"]=aa
}


AA=0
aa=0
Aa=0
i=0
j=0

## populate in df the count of AA,Aa, and aa for pd
for (i in 1:nrow(pd_map)) {
  
  
  Aa=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  AA=sum(pd_map[i,4]==pd_pre[,2*i+1] & pd_map[i,4]==pd_pre[,2*i+2]) 
  aa=sum(pd_map[i,5]==pd_pre[,2*i+1] & pd_map[i,5]==pd_pre[,2*i+2]) 
  
  df[i,"AA_pd"]=AA
  df[i,"Aa_pd"]=Aa
  df[i,"aa_pd"]=aa
  
}
#populate in df the chromosome and locus

df[,1]=cc_map[,1]
df[,2]=cc_map[,3]


##Store the dataframe df for chromosome XY
dfXY=df
############################################################


###Row bind the dataframes for a cpmplete dataset for all chromosomes for all locii
df_combined= rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22,dfX,dfY,dfXY)

############################
############################
##first define the function to apply
Chsq <- function(x){
  ## input is a row of your data
  ## creating a table from each row
  if (x[1]!=0 & x[2]!=0 & x[3]==0)  y <- matrix(x[c(1,2,4,5)],byrow =TRUE,nrow=2)
  if (x[1]!=0 & x[2]!=0 & x[3]!=0)  y <- matrix(x[c(1,2,3,4,5,6)],byrow =TRUE,nrow=2)
  if (x[1]==0 & x[2]!=0 & x[3]!=0)  y <- matrix(x[c(2,3,5,6)],byrow =TRUE,nrow=2)
  if (x[1]==0 & x[2]==0 & x[3]!=0)  y <- matrix(x[c(3,6)],byrow =TRUE,nrow=2)
  if (x[1]==0 & x[2]==0 & x[3]==0)  return (1) ##for all control zero count 
  if (x[1]!=0 & x[2]==0 & x[3]!=0)  y <- matrix(x[c(1,3,4,6)],byrow =TRUE,nrow=2)
  if (x[1]==0 & x[2]!=0 & x[3]==0)  y <- matrix(x[c(2,5)],byrow =TRUE,nrow=2)
  if (x[1]!=0 & x[2]==0 & x[3]==0)  y <- matrix(x[c(1,4)],byrow =TRUE,nrow=2)
  
  ### this will return the p value
  return(chisq.test(y)$p.value)
}
data =df_combined[,3:8]
P_Values <- as.vector(apply(data,1,Chsq))   ###return p values
df_combined = cbind(df_combined,P_Values) ####Populate p values in the df 


######################################


#sort the combined data frame for all chromoaomes based on lowest p-value value to lowest
sorted_df_comb_by_p_values= df_combined[with(df_combined, order(P_Values)), ]







###Q1.
#Top ten SNP Lociassociated with Parkinson's Disease, ordered by chi-sq test p-values
sorted_df_comb_by_p_values[1:10,]

###Q3
######################### #SNP Locus rs3741411###################################


sorted_df_comb_by_p_values[sorted_df_comb_by_p_values$locus == "rs3741411", ]
y=matrix(c(184,13,161,9),byrow =TRUE,nrow=2)
y
chisq.test(y)$p.value
###Q4 Number of SNP loci which have a p-value of <0.05
nrow(sorted_df_comb_by_p_values[sorted_df_comb_by_p_values$P_Values<0.05,])
