###################################################
############ 20240123 update by GR ################
############ bulk RNA-seq analysis ################
#################### for DEG ######################
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
############# heatmap 구현#############
#######################################

# heatmap에 표현할 유전자를 지정(보고싶은 gene 이 있으면 그걸 지정해주면 됨 
# (**** 주의!!!! ****  input시킨 excel raw file에 있는 gene 이름 그대로 입력해야함 대소문자도 정확히!)
# (이 과정은 필요에 따라 진행 만약 모든 gene을 다 사용할 경우 이 과정 필요 없음.)
selected.genes = c("CRP","COL5A3","KCNB1","SLC16A14","GOLM1","CSTA","SAA1","ACE2","FABP4","ITLN1","LDHC","ARHGEF37","CYP26A1","MEP1B","OLFM2","CCRN4L","ARNTL","LPL","SELE","ID1","GNAO1","MLIP","LINC00261","UGT1A8","COL4A5","SDCBP2","FMO1","SEPT4","ZCCHC16","ABCB4","NPAS2","PRKCE","UPP2","TMTC1","CLDN2","UBXN10","FLG","FABP5","B3GAT1","BST1","COL28A1","DDX11","EYA2","FKBP5","GREM1","GSTM3","NR1D1","NR1D2","NUDT13","PCOLCE","PLEKHA6","SCAMP5","SLC2A4","SNHG5","TNIK","TUBA4A","TYW1B")
selected.genes.df = data.norm.table[ selected.genes, ]

# normalization(expression수치 범위 너무 넓어서) + 전체 합 0인 행 제거하기
norm_selected.genes.df <- log10(selected.genes.df + 1)
selected.genes.df <- norm_selected.genes.df[which(rowSums(norm_selected.genes.df) != 0),]

# heatmap color 지정(바꾸려면 컬러차트 검색해서 추가, 제외 하면 됨.)
hmcol <- colorRampPalette(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC"))(n = 10)
# heatmap 그리기
library(pheatmap)
pheatmap(selected.genes.df, scale = "row")
