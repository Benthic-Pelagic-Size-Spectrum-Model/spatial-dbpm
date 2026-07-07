ad <- read.csv("all_dint.csv")   # lme, mask_id, int200, int_bw(old area-mean-of-log), dint
out <- data.frame(lme=integer(), dint=numeric())
for (i in seq_len(nrow(ad))) {
  L <- ad$lme[i]; f <- sprintf("int_lme%d.csv", L)
  if (!file.exists(f)) { out <- rbind(out, data.frame(lme=L, dint=ad$dint[i])); next }
  d <- read.csv(f); names(d) <- tolower(names(d)); d <- d[is.finite(d$int_bw), ]
  if (nrow(d) < 3) { out <- rbind(out, data.frame(lme=L, dint=ad$dint[i])); next }
  w <- cos(d$lat*pi/180); b <- 10^d$int_bw
  int_bw_hw <- log10(sum(w*b)/sum(w))          # arithmetic biomass mean = average-then-fit intercept
  dint_new  <- int_bw_hw - ad$int200[i]        # total offset from parquet baseline (~int200)
  out <- rbind(out, data.frame(lme=L, dint=round(dint_new, 3)))
}
write.csv(out, "lme_dint_hw.csv", row.names=FALSE)
# report the change vs old dint
ad$dint_new <- out$dint[match(ad$lme, out$lme)]; ad$gain <- round(ad$dint_new - ad$dint, 2)
cat("median gain over old dint:", median(ad$gain,na.rm=TRUE), " | weak regions:\n")
w <- ad[ad$lme %in% c(171,131,177,134,161,62,16,55,63,66,23), c("lme","dint","dint_new","gain")]
print(w, row.names=FALSE)
cat("sample good LMEs:\n"); print(ad[ad$lme %in% c(1,3,5,22), c("lme","dint","dint_new","gain")], row.names=FALSE)
