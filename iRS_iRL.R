### DATA

#Communities for each treatment
Co.dat <- c(6.409675, 6.194636, 6.452469, 6.072689, 6.275555, 6.263966, 6.103991, 6.440885, 6.637393, 6.579923)
Cx.dat <- c(6.440099, 5.774942, 6.129712, 5.993228, 6.313750)        
Po.dat <- c(2.815681, 5.858012, 4.413149, 5.412506, 4.225717, 3.606892, 2.550436, 2.820317, 2.054459)
Px.dat <- c(6.258275, 6.573278, 6.384016, 5.890846, 6.250952, 6.729738, 6.419818, 6.331941, 5.901057)

#Data could be in a dataframe
shannon<- matrix(c(Co.dat,Po.dat,Cx.dat,Px.dat),nrow=1) ; 
my.dat<-shannon
row.names(my.dat) <- "shannon"  # Other variables could be included in the matrix
my.dat.fra<-rbind(my.dat,my.dat,my.dat); #Here we included three times the same var
row.names(my.dat.fra) <- c("Var1","Var2","Var3")
#we need  a factor
my.factor<-c(rep("Co",length(Co.dat)),rep("Po",length(Po.dat)),rep("Cx",length(Cx.dat)),rep("Px",length(Px.dat)))
my.factor<-factor(my.factor, levels=c("Co","Po","Cx","Px"))


## Load a  the function first
# -for one variriable
fun.RS.RL.BOOT.2(shannon, "Shannon") # the output is in clipboard
fun.RS.RL.BOOT.2(my.dat.fra[3,], rownames(my.dat.fra)[3])

# -for a dataframe with variables in rows
out<- fun.RS.RL.BOOT.2(my.dat.fra[1,], rownames(my.dat.fra)[1]) # we need the out structure
for (i in 2:nrow(my.dat.fra)) { 
   out[i,]<-fun.RS.RL.BOOT.2(my.dat.fra[i,], rownames(my.dat.fra)[i])}
out
# -write output to a file
write.table(out, "RS_RL_out.txt", quote=F, sep="\t", row.names=T, col.names=T)


# - A plot from one variable

fun.RS.RL.BOOT.2(shannon, "Shannon")
fun.RL.plot(shannon, "Shannon")



## Functions
fun.RS.RL.BOOT.2<-function(y,tit){
  y<-c(t(y));
  nboot<-100; nsamples=10   # The number of replications and samples could be modified 
  my.mod<-data.frame(matrix(c(rep(1,25)),nrow=1,byrow=T)); rownames(my.mod)<-tit
  colnames(my.mod)<-c("Co","Po","Cx","Px","iRS","iRL.OW","iRL","Co-Cx",
                    "p(Co=Po).ml","p(Co=Cx).ml","p(Co=Px).ml",   # ml check-points  9-11
                    "m(iRS)","sd(iRS)","p(iRS=1)",               # iRS boot        12-14
                    "m(iRL.OW)","sd(iRL.OW)","p(iRL.OW=1)",      # iRL.OW boot     15-17
                    "m(iRL)","sd(iRL)","p(iRL=1)",               # iRL boot        18-20
                    "m(Co-Cx)","sd(Co-Cx)","p((Co-Cx)=0)",       # Co-Cx boot      21-23
                    "Clasif.ml","Clasif.boot")                   #                 24-25
  my.mod[1:4]<-tapply (y, my.factor,mean)
  Co<-y[which(my.factor=="Co")]; Po<-y[which(my.factor=="Po")]
  Cx<-y[which(my.factor=="Cx")]; Px<-y[which(my.factor=="Px")]
  my.mod[5:7]<-fun.RS.RL(Co,Cx,Po,Px); my.mod[8]<-mean(Co)-mean(Cx)
  # ml check-points
      CoPo<-which(my.factor %in% c("Co","Po"))
      my.mod[9]<- anova(lm(y[CoPo]~ my.factor[CoPo]))$"Pr(>F)"[1];  if(my.mod[9]=="NaN"){my.mod[9]<-1}
      CoCx<-which(my.factor %in% c("Co","Cx"))
      my.mod[10]<- anova(lm(y[CoCx]~ my.factor[CoCx]))$"Pr(>F)"[1];  if(my.mod[10]=="NaN"){my.mod[10]<-1}
      CoPx<-which(my.factor %in% c("Co","Px"))
      my.mod[11]<- anova(lm(y[CoPx]~ my.factor[CoPx]))$"Pr(>F)"[1];  if(my.mod[11]=="NaN"){my.mod[11]<-1}
      my.mod[24]<-fun.clasif(as.numeric(my.mod[c(9:11,7)]),as.numeric(my.mod[2])-as.numeric(my.mod[1]))
  # bootstrap   
    my.boot<-matrix(c(rep(1,4*nboot)),nrow=nboot,byrow=T)
    for (ii in 1:nboot){
      Co.B <-mean(sample(Co,10,r=T))
      Cx.B <-mean(sample(Cx,10,r=T))
      Po.B <-mean(sample(Po,10,r=T))
      Px.B <-mean(sample(Px,10,r=T)) 
      my.boot[ii,1:3] <- fun.RS.RL(Co.B,Cx.B,Po.B,Px.B)
      my.boot[ii,4]<-1+mean(Co.B)-mean(Cx.B)
    }
  for(ii in 1:4){
    meB<-mean(my.boot[,ii]); sdB<-sd(my.boot[,ii])
    cola=0; if(meB >=1) {cola=1}
    PrB<-pnorm(1,meB,sdB,lower.tail = cola)
    my.mod[(ii*3+9):(ii*3+11)]<-c(meB,sdB,PrB)
  }
  my.mod[21]<-my.mod[21]-1
  my.mod[25]<-fun.clasif(as.numeric(my.mod[c(14,23,20,18)]),as.numeric(my.mod[2])-as.numeric(my.mod[1]))   ### alternativa cambiar 18 por 7
  write.table(my.mod, "clipboard", quote=F, sep="\t", row.names=T, col.names=T)
  my.mod
}

fun.clasif<-function(cc, rr){  #  (P.RS, P.Cx-Co, P.RL, iRL, Po-Co)
  pert<-ifelse (rr > 0 ,"+",  "-") ; myclasif= "RS"
  if (cc[1] >= .05) {myclasif= "RS";  pert =""}
  else {
    if (cc[2] < 0.01) {myclasif= "TV"}
    else{
      if (cc[3] >= .05) {myclasif = "RL"}
      else{
        if ((cc[4] > 0) & (cc[4] < 1)) {myclasif = "RC"}
        if (cc[4] <= 0)  {myclasif = "CT"}
        if (cc[4] > 1)   {myclasif = "RB"}
      } } }
  myclasif<-paste(myclasif, pert, sep="")
  myclasif
}



fun.RS.RL<-function(Co,Cx,Po,Px){ 
  mCo <- mean(Co); mCx <- mean(Cx); mPo <- mean(Po); mPx <- mean(Px)
  Do <-abs(mCo-mPo);  Dx <-abs(mCx-mPx) 
  if ((mCo == 0 & mPo == 0) | (mCo ==  mPo)) {RS=1; RL.OW=1; RL =1}
  else {RS <- 1-2*Do/(mCo+Do); RL.OW <- 2*Do/(Do+Dx)-1 ; RL <- (mPx-mPo)/(mCo-mPo)}
  c(RS, RL.OW, RL)
}

fun.RL.plot<-function(y,tit){
  y<-c(t(y));
  boxplot(y~my.factor,col=c("green", "red", "lightgreen", "orange"), main=tit, xlab="", ylab="")
 # stripchart(my.dat[i,]~my.factor, vertical = TRUE, method = "jitter", add = TRUE, pch = 16)
  }


