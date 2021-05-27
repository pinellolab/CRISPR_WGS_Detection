library(RColorBrewer)
library('data.table')


cols = c(brewer.pal(9,'Set1'),'black')


plotGuideseq = function(guideName)
{
	d = fread(paste0(guideName,'.guideseq.txt'))
	dd = d[,c(1,2,3,4,13,19,22,23,24,32)]
	colnames(dd) = c('g.chr','g.start','g.end','g.guideseq','g.guideseq_reads','g.mismatches','w.chr','w.start','w.end','w.pct')
#	head(dd[dd$w.chr != ".",])

	dd$wgsVal = as.numeric(dd$w.pct)
	dd[is.na(dd$wgsVal),]$wgsVal = 0
	dd$col = dd$g.mismatches + 1
	dd[dd$col > 9,]$col = 10
	plot(dd$wgsVal,dd$g.guideseq_reads,xlab="WGS Indel Ratio",ylab="GUIDE-seq Reads",col=cols[dd$col],pch=20,main=paste(guideName,'\nN =',nrow(dd),'GUIDE-seq hits'))
	col_labs = as.numeric(names(table(dd$g.mismatches)))+1
	col_names = col_labs - 1
	legend('top',fill=cols[col_labs],legend=col_names,title='#Mismatches',ncol=2)
}

plotWGS = function(guideName)
{
	d = fread(paste0(guideName,'.wgsHits.txt'))
	d = d[d$V7 > 100 & d$V10 > 100,]
	dd = d[,c(1,2,3,11,12,14,13,15,24,30)]
	colnames(dd) = c('w.chr','w.start','w.end','w.pct','g.chr','g.start','g.end','g.guideseq','g.guideseq_reads','g.mismatches')
#	head(dd[dd$w.chr != ".",])
	dd = dd[order(dd$w.pct,decreasing=T)[1:min(nrow(dd),2500)],]

	dd$guideseqVal = as.numeric(dd$g.guideseq_reads)
	dd[is.na(dd$guideseqVal),]$guideseqVal = 0
	dd$col = as.numeric(dd$g.mismatches)
	dd[is.na(dd$col),]$col = 9
	dd$col = dd$col+1
	plot(dd$w.pct,dd$guideseqVal,xlab="WGS Indel Ratio",ylab="GUIDE-seq Reads",col=cols[dd$col],pch=20,main=paste(guideName,'\nShowing the top',nrow(dd),'WGS hits'))
	col_labs = as.numeric(setdiff(names(table(dd$g.mismatches)),"."))+1
	col_names = col_labs - 1
	legend('top',fill=cols[c(col_labs,10)],legend=c(col_names,'WGS only'),title='#Mismatches',ncol=2)
}

for (guide in c('Test'))
{
	plotGuideseq(guide)
	plotWGS(guide)
}

