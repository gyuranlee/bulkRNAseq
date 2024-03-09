###################################################
############ 20240123 update by GR ################
############ bulk RNA-seq analysis ###############
###################################################

###### 필요한 package 설치 ######
# R에서 필요한 아래의 2가지 package를 install 해야함.
# 1. BiocManager
# 2. Deseq2
# 위의 2개 package를 인스톨하기 위해 아래의 명령어 실행.(#이 명령어 앞에 붙어있으면 실행되지 않아요! 만약 사용자의 R에 아래의 package
# 가 설치되어있지 않은 상태라면 아래 명령어 다섯줄 앞에 붙은 # 표시를 지우고 실행해 주세요:)  )
################################

#install.packages("BiocManager")
#BiocManager::install("DESeq2", force=TRUE)
#install.packages("ggplot2")
#install.packages('reshape2')
#sessionInfo()


# R로 방금 설치한 package를 불러오는 작업(library함수로 package 불러오는 작업은 매 실행시마다 해야함)
library("DESeq2")
library("ggplot2")
library("reshape2")

# 현재 working 디렉토리를 설정
# 사용자가 원하는 자신의 컴퓨터 디렉토리를 지정하시면 됨(raw count table이 저장되어있는 경로로 ''안에 작성하면 됨)
# 만약 현재 디렉토가 어디인지 확인하고 싶으시면 getwd() 명령어를 사용하세요!
setwd('~/DIRECTORY')

# raw count table을 불러옵니다.
raw.count = read.csv(file="IFT46kD_in_IMCD/IMCD_count/mouse_cell_line_raw.csv", header=TRUE)

# 불러온 테이블 확인
dim(raw.count) # row, column갯수에 대한 정보 확인 가능
head(raw.count)

# DESeq2가 원하는 형식으로 테이블을 만들어 줍니다.
# 테이블의 Row name 을 유전자로 만들어주고, 불필요한 column은 제외시켜줍니다.
rownames(raw.count) = make.names(raw.count[,1], unique=TRUE)
raw.count = raw.count[,-c(0:1)]

# 수정한 테이블 확인
dim(raw.count) 
head(raw.count)

# 이제 DESeq2에 input으로 들어갈 table에 대한 정보를 담고 있는 sample.info 데이터 프레임을 만들어줌.
# (어떤열이 control인지 어떤게 질병인ㄷ지에 대해 정리할 dataframe을 하나 더 만들어야함.)
# *** 이때 주의할 점! 위의 count table과 column명을 동일하게 만들어야 함 ***
sample.info = data.frame(row.names=c("CTRL_1","CTRL_2","CTRL_3","IFT46_1","IFT46_2","IFT46_3"),
                         condition= factor(c("Control","Control","Control","Sample","Sample","Sample")))

# 방금 만든 sample.info 테이블 확인
dim(sample.info)
head(sample.info)

##########################################################
###### DESeq package에 들어갈 데이터 모두 준비 완료 ######
##########################################################

# sample.info + raw count --> 'DESeqDataSetFromMatrix' function을 통해 dataset으로 변환
data = DESeqDataSetFromMatrix(countData = round(raw.count), colData = sample.info ,design =~ condition)

# 변환된 dataset에 대한 정보를 살펴봅니다.
data

# Normalization
# (컴퓨터 환경에 따라 다르긴 하지만, 약 15-20초 소요.)
data = DESeq(data) # normalization function

# Normalization 되어진 table 불러옴.
data.norm.table = counts(data, normalized=T)

# normalization 하기 전과 후의 table을 비교
head(raw.count) # normalization 전
dim(raw.count) # normalization 전

head(data.norm.table) # normalization 후
dim(data.norm.table) # normalization 후


#######################################
#############   for DEG  ##############
#######################################

# DEG 찾기
# (results 라는 function 사용 :: input file - normalization 한 gene expression data)
res <- results(data)

# DEG table + expression table 합치기
res.exp <- cbind(res, data.norm.table)
head(res.exp)

# p-value 순 정렬 (NA값인 유전자 제외)
res.exp <- res.exp[order(res.exp$pvalue),]
res.filt <- res.exp[!is.na(res.exp$padj),]

# Q-value와 log2FC 범위 설정 ( 필요에 따라 수치 변경하시면 됩니다!)
Qvaule = 0.05
log2FC = 1

#p-value 0.05보다 낮고, log2foldchange 1보다 낮은거 동시에 필터링 (&쓰면 교집합의미)
up.regulated.genes <- res.filt[ which( (res.filt$padj <= Qvaule) & (res.filt$log2FoldChange > log2FC) ) ,]
down.regulated.genes <- res.filt[ which( (res.filt$padj <= Qvaule) & (res.filt$log2FoldChange < -log2FC) ) ,]

dim(up.regulated.genes)
dim(down.regulated.genes)

## 필요에 따라 아래 gene list 파일로 저장하면 됨 (**앞서 지정해둔 경로에 파일로 저장됨**) ##
# up-regulated gene 파일로 저장 
write.table(up.regulated.genes, file ="Up_regulated_genes_result.txt", row.names = T, col.names = T, sep ="\t", quote = F)
# down-regulated gene 파일로 저장 
write.table(down.regulated.genes, file ="Down_regulated_genes_result.txt", row.names = T, col.names = T, sep ="\t", quote = F)
# 전체 gene 파일로 저장
write.table(res.exp, file ="All_genes_result.txt", row.names = T, col.names = T, sep ="\t", quote = F)


##### visualization #####
# (필요에 따라 진행하면 됨 : DEG 결과를 visualization 해서 확인하고 싶은 경우)

# volcano plot
ggplot(data.frame(res.filt), aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(color='#37465d', alpha=0.5) + #알파값은 투명도
  geom_hline(yintercept=-log10(Qvaule), linetype='dashed', color='#f2c1b9') +
  geom_vline(xintercept=c(-log2FC,log2FC), linetype='dashed', color='#f2c1b9') +
  geom_text(aes(label=ifelse(padj < res.filt[11, "padj"], rownames(res.filt)[1:10],'')), hjust=0.5, vjust=-0.5) +
  #top 10개만 이름을 가져오겠다. (hjust는 글씨의 위치 )
  theme_classic()

