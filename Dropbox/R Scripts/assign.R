design<-model.matrix(~0+factor(c(1,1,1,1,1,1,2,2,2,2,2)))
colnames(design) <- c("NormDiet", "KetoDiet")
contrast.matrix <- makeContrasts(NormDiet-KetoDiet,levels=design)
contrasts.list = c("NormDiet vs KetoDiet")


fit.nobg<-lmFit(eset.nobg.constant,design)
fit.nobg.contrast<-contrasts.fit(fit.nobg,contrast.matrix)
fit.nobg.contrast<-eBayes(fit.nobg.contrast)
topTable(fit.nobg.contrast, coef=1, adjust = "BH", n=40)

fit.bg<-lmFit(eset.bg.constant,design)
fit.bg.contrast<-contrasts.fit(fit.bg,contrast.matrix)
fit.bg.contrast<-eBayes(fit.bg.contrast)
topTable(fit.bg.contrast, coef=1, adjust = "fdr", n=20)

fit.nobg.contrast[fit.nobg.contrast$gene == "1370111_at",]

Table<-topTable(fit.nobg.contrast, coef=1, adjust = "BH", n=nrow(fit.nobg.contrast))

Table[Table$ID == "1370111_at",]
