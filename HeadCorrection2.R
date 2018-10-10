require(tidyverse)
install.packages("pracma")
require(pracma)
EMAraw<-read_tsv(file.choose())

#TTraw (Tongue Tip) has coordinates (Tx_4, Ty_4, Tz_4)
#but should be redefined from origin at UI.
EMAraw$TTrawx = EMAraw$Tx_4-EMAraw$Tx_6
EMAraw$TTrawy = EMAraw$Ty_4-EMAraw$Ty_6
EMAraw$TTrawz = EMAraw$Tz_4-EMAraw$Tz_6
#TDraw (Tongue Dorsum) has coordinates (Tx_3, Ty_3, Tz_3)
#but should be redefined from origin at UI.
EMAraw$TDrawx = EMAraw$Tx_3-EMAraw$Tx_6
EMAraw$TDrawy = EMAraw$Ty_3-EMAraw$Ty_6
EMAraw$TDrawz = EMAraw$Tz_3-EMAraw$Tz_6
#LIraw (Lower Incisors) has coordinates (Tx_5, Ty_5, Tz_5)
#but should be redefined from origin at UI.
EMAraw$LIrawx = EMAraw$Tx_5-EMAraw$Tx_6
EMAraw$LIrawy = EMAraw$Ty_5-EMAraw$Ty_6
EMAraw$LIrawz = EMAraw$Tz_5-EMAraw$Tz_6

EMAnorm<-select(EMAraw,
                time_s = `Wav Time`,
                #A (Bridge of Nose) has coordinates (Tx, Ty, Tz)
                Ax = Tx,
                Ay = Ty,
                Az = Tz,
                #B (Left Mastoid) has coordinates (Tx_1, Ty_1, Tz_1)
                Bx = Tx_1,
                By = Ty_1,
                Bz = Tz_1,
                #C (Right Mastoid) has coordinates (Tx_2, Ty_2, Tz_2)
                Cx = Tx_2,
                Cy = Ty_2,
                Cz = Tz_2,
                #UI (Upper Incisors - Maxilla) has coordinates (Tx_6, Ty_6, Tz_6)
                UIx = Tx_6,
                UIy = Ty_6,
                UIz = Tz_6,
                #TTraw (Tongue Tip) has coordinates (Tx_4, Ty_4, Tz_4)
                #but should be redefined from origin at UI.
                TTrawx = TTrawx,
                TTrawy = TTrawy,
                TTrawz = TTrawz,
                #TDraw (Tongue Dorsum) has coordinates (Tx_3, Ty_3, Tz_3)
                #but should be redefined from origin at UI.
                TDrawx = TDrawx,
                TDrawy = TDrawy,
                TDrawz = TDrawz,
                #LIraw (Lower Incisors) has coordinates (Tx_5, Ty_5, Tz_5)
                #but should be redefined from origin at UI.
                LIrawx = LIrawx,
                LIrawy = LIrawy,
                LIrawz = LIrawz)

#Matrices of coordinates. One row = one frame
A<-cbind(EMAnorm$Ax,EMAnorm$Ay,EMAnorm$Az)
B<-cbind(EMAnorm$Bx,EMAnorm$By,EMAnorm$Bz)
C<-cbind(EMAnorm$Cx,EMAnorm$Cy,EMAnorm$Cz)
UI<-cbind(EMAnorm$UIx,EMAnorm$UIy,EMAnorm$UIz)
TTraw<-cbind(EMAnorm$TTrawx,EMAnorm$TTrawy,EMAnorm$TTrawz)
TDraw<-cbind(EMAnorm$TDrawx,EMAnorm$TDrawy,EMAnorm$TDrawz)
LIraw<-cbind(EMAnorm$LIrawx,EMAnorm$LIrawy,EMAnorm$LIrawz)
#M is the projection of A to line BC with coordinates (Mx, My, Mz)
#Matrices of vectors. One row = one frame
vBC<-cbind(C-B)
vBA<-cbind(A-B)
#Transpose vector matrices because dot() function operates on matrix columns!
w<-dot(t(vBA),t(vBC))/dot(t(vBC),t(vBC))
M = B + w * vBC

#Vectors vMA (posterior-anterior) and vMB (lateral) should be perpendicular to one another
vMA<-cbind(A-M)
vMB<-cbind(B-M)

#How many frames have vMA and vBC perfectly perpendicular?
cat("The Posterior-Anterior and Lateral axes are perfectly perpendicular in",
    sum((dot(t(vMA),t(vBC)))==0,na.rm=TRUE),
    "out of",
    length((dot(t(vMA),t(vBC)))),
    "frames.")

####A handful of these should be exactly 0.
###The rest will be slightly off due to floating point stuff.
###Think it's not a huge issue

#vNorm is the cross product of vMA and vMB and perpendicular to both
EMAnorm$vNormi = (EMAnorm$vMAj*EMAnorm$vMBk - EMAnorm$vMAk*EMAnorm$vMBj)
EMAnorm$vNormj = (EMAnorm$vMAk*EMAnorm$vMBi - EMAnorm$vMAi*EMAnorm$vMBk)
EMAnorm$vNormk = (EMAnorm$vMAi*EMAnorm$vMBj - EMAnorm$vMAj*EMAnorm$vMBi)
vNorm<-cbind(EMAnorm$vNormi,EMAnorm$vNormj,EMAnorm$vNormk)
vNorm<-c(vNormi,vNormj,vNormk)

#Normalise basis vectors to unit vectors. MAKE SURE TO USE rowSums AND NOT sum!
vPostAnt = vMA / (sqrt(rowSums(vMA*vMA,na.rm=TRUE)))
vLateral = vMB / (sqrt(rowSums(vMB*vMB,na.rm=TRUE)))
vSupInf = vNorm / (sqrt(rowSums(vNorm*vNorm,na.rm=TRUE)))

#Redefine TTraw according to new axes, this is TT
EMAnorm$TTx<-dot(t(TTraw),t(vPostAnt))
EMAnorm$TTy<-(dot(t(TTraw),t(vSupInf)))*(-1)
EMAnorm$TTz<-dot(t(TTraw),t(vLateral))
TT<-c(TTx,TTy,TTz)
#Redefine TDraw according to new axes, this is TD
EMAnorm$TDx<-dot(t(TDraw),t(vPostAnt))
EMAnorm$TDy<-(dot(t(TDraw),t(vLateral)))*(-1)
EMAnorm$TDz<-dot(t(TDraw),t(vSupInf))
TD<-c(TDx,TDy,TDz)
#Redefine LIraw according to new axes, this is LI
EMAnorm$LIx<-dot(t(LIraw),t(vPostAnt))
EMAnorm$LIy<-(dot(t(LIraw),(vLateral)))*(-1)
EMAnorm$LIz<-dot(t(LIraw),t(vSupInf))
LI<-c(LIx,LIy,LIz)

#Plot tongue tip height with red line at upper incisor sensor height
TTHeight<-ggplot(EMAnorm, aes(time_s, TTy, alpha=TTx))
TTHeight + geom_line(size=1) + geom_hline(aes(yintercept=0), colour="red") + theme_bw() +
  scale_x_continuous(limits = c(3.9, 4.37)) + scale_y_continuous(limits = c(-15, 15)) +
  ylab("") + xlab("") + theme(legend.position="none")

