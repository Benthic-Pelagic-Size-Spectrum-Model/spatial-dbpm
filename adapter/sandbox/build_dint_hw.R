suppressMessages(library(arrow))
ad<-read.csv("all_dint.csv"); out<-data.frame(lme=integer(),dint=numeric())
for(i in seq_len(nrow(ad))){ L<-ad$lme[i]; f<-sprintf("int_lme%d.csv",L)
  pf<-Sys.glob(sprintf("DBPM_dev_uv/dbpm_inputs/dbpm_clim-fish-inputs_fao_lme-%d_*.parquet",L))[1]
  if(!file.exists(f)||is.na(pf)){ out<-rbind(out,data.frame(lme=L,dint=ad$dint[i])); next }
  d<-read.csv(f); names(d)<-tolower(names(d)); d<-d[is.finite(d$int_bw),]
  if(nrow(d)<3){ out<-rbind(out,data.frame(lme=L,dint=ad$dint[i])); next }
  w<-cos(d$lat*pi/180); b<-10^d$int_bw; ibw<-log10(sum(w*b)/sum(w))
  pint<-mean(read_parquet(pf,col_select="intercept")$intercept,na.rm=TRUE)
  out<-rbind(out,data.frame(lme=L,dint=round(ibw-pint,3))) }
write.csv(out,"lme_dint_hw.csv",row.names=FALSE)
cat(sprintf("wrote lme_dint_hw.csv: %d regions, %d NA remaining\n", nrow(out), sum(is.na(out$dint))))
cat("LME 32 dint now:", out$dint[out$lme==32], "\n")
