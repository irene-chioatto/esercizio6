library(fImport)
library(GAS)
library(backtest)
library(zoo)
library(FinTS)
library(forecast)
library(ggplot2)
library(fGarch)
library(lmtest)
library(rugarch)


data<-fredSeries("DEXUSEU",from="1999-01-04",to="2020-05-01")
data<-as.zoo(data)
dates<-as.Date(index(data))
dex<-as.numeric(data)
n<-length(dex)
dexuseu<- data.frame(
  day = dates[-1],
  value = dex[-1], 
  lvalue = log(dex)[-1],
  lret = log(dex[2:n]) - log(dex[1:(n-1)]),
  lret.sq = (log(dex[2:n]) - log(dex[1:(n-1)]))^2,
  lret.abs = abs(log(dex[2:n]) - log(dex[1:(n-1)])))  

attach(dexuseu)

##Punto 1. Media condizionata
ggplot(data=dexuseu, aes(x=day, y=value)) + geom_line(color="deeppink2") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels="%d %b %Y") + labs(title = "DEXUSEU", subtitle="Serie da 04-01-1999 a 01-05-2020", x = "time", y = "index")
#Dal grafico si nota una evidente non stazionariet? della serie considerata. Presenta un trend crescente dal 2002 al 2008, mentre dopo 
#la crisi del 2008 il trend diventa pressocch? decrescente.
require(gridExtra)
plot1<-ggAcf(value, type = c("correlation"),lag.max=100) +labs(title = "ACF", x = "lag", y = "ACF")
plot2<-ggPacf(value,lag.max=100) + labs(title = "PACF", x = "lag", y = "PACF")
grid.arrange(plot1, plot2, nrow=2)
#L'ACF conferma la presenza di dipendenza temporale e la necessit? di differenziare la serie.
ggplot(data=dexuseu, aes(x=day, y=lret)) + geom_line(color="deeppink2") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels="%d %b %Y") + labs(title = "Log-returns di DEXUSEU", subtitle="Serie da 04-01-1999 a 01-05-2020", x = "time", y = "index")
# statistiche pi? importanti dei log rendimenti
summary(lret)
round(basicStats(lret), 4)
#possiamo notare che il valore atteso del rendimento ? zero.
#la serie risulta essere leggermente asimmetrica a destra, inoltre le code della distribuzione sono poco pi? pesanti della Normale.
#possiamo verificarlo anche guardando al grafico quantile-quantile:
qqnormPlot(lret)
#Analizziamo il correlogramma:
require(gridExtra)
plot1<-ggAcf(lret, type = c("correlation"),lag.max=200) +labs(title = "ACF dei log-returns", x = "lag", y = "ACF")
plot2<-ggPacf(lret,lag.max=100) + labs(title = "PACF dei log-returns", x = "lag", y = "PACF")
grid.arrange(plot1, plot2, nrow=2)
#Si notano alcuni ritardi significativi che si protraggono a lungo nel tempo e che non sembrano seguire un particolare
#pattern che ci possa far pensare ad una stagionalit?.
#Quest'effetto ? dovuto con tutta probabilit? all'informazione non modellata circa l'eteroschedasticit?, nonostante 
#stiamo considerando la serie al grado 1.




#Punto 2. Varianza condizionata
#Dal grafico dei log-returns del DEXUSEU ? evidente il fenomeno del volatility clustering, segno della presenza di eteroschedasticit?.
#valutiamo quindi i log-rendimenti al quadrato e i in valore assoluto:
ggplot(data=dexuseu, aes(x=day, y=lret.sq)) + geom_line(color="deeppink2") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels="%d %b %Y") + labs(title = "Log-returns al quadrato", subtitle="Serie da 04-01-1971 a 01-05-2020", x = "time", y = "index")

plot1<-ggAcf(lret.sq, type = c("correlation"),lag.max=100) +labs(title = "ACF dei log-returns al quadrato", x = "lag", y = "ACF")
plot2<-ggPacf(lret.sq,lag.max=100) + labs(title = "PACF dei log-returns al quadrato", x = "lag", y = "PACF")
grid.arrange(plot1, plot2, nrow=2)
#si nota una forte dipendenza nella volatilit? dei rendimenti
AutocorTest(lret.sq)
ArchTest(lret.sq)

ggplot(data=dexuseu, aes(x=day, y=lret.abs)) + geom_line(color="deeppink2") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels="%d %b %Y") + labs(title = "Log-returns in valore assoluto", subtitle="Serie da 04-01-1971 a 01-05-2020", x = "time", y = "index")
plot1<-ggAcf(lret.abs, type = c("correlation"),lag.max=100) +labs(title = "ACF dei log-returns in valore assoluto", x = "lag", y = "ACF")
plot2<-ggPacf(lret.abs,lag.max=100) + labs(title = "PACF dei log-returns in valore assoluto", x = "lag", y = "PACF")
grid.arrange(plot1, plot2, nrow=2)
#anche valutando il valore assoluto, il correlogramma mostra dipendenza.
AutocorTest(lret.abs) 
ArchTest(lret.abs)
#? evidente la presenza di effetti ARCH sia sulla serie dei log-rendimenti al quadrato che in valore assoluto.



#Punto 3.
skewness(lret)
kurtosis(lret)
#come accennato al punto precedente, la serie presenta una leggera asimmetria a destra.
#l'eccesso di curtosi positivo inoltre indica una distribuzione leptocurtica.
#possiamo verificare la deviazione dalla Normale guardando al grafico quantile-quantile e ad alcuni test di normalita':
qqnormPlot(lret)
#i dati si discostano dalla bisettrice in particolare sulla coda destra.
ksnormTest(lret)
jarqueberaTest(lret)
dagoTest(lret)
#proviamo ad adattare un modello sulla distribuzione:
(fit.norm<-nFit(lret))
(fit.std<-stdFit(lret)) #t di Student
(fit.ged<-gedFit(lret)) #GED
(fit.sstd<-sstdFit(lret)) #t-Student asimmetrica
(fit.sged<-sgedFit(lret)) #GED asimmetrica
(fit.nig<-nigFit(lret,trace=F)) #NIG
N<-length(lret)
#calcoliamo i valori del criterio di Akaike per ogni distribuzione.
(aic.norm<-2*fit.norm@fit$minimum + 2*2/N)
(aic.std<-2*fit.std$objective + 2*3*N)
(aic.ged<-2*fit.ged$objective + 2*3*N)
(aic.sstd<-2*fit.sstd$minimum + 2*4*N)
(aic.sged<-2*fit.sged$objective + 2*4*N)
(aic.nig<-2*fit.nig@fit$objective + 2*4*N)
which.min(c(aic.norm,aic.std,aic.ged,aic.sstd,aic.sged,aic.nig))
#viene preferito il secondo modello inserito, ossia una distribuzione student-t
#densita' empirica a confronto con quella teorica:
ggplot(data=dexuseu, aes(x=lret)) + geom_histogram(aes(y=..density..),bins = 100, col="white", fill="steelblue") +
  geom_density(col=3,lwd=1) + theme(axis.text.x=element_text(angle=0, hjust=1, size=16), axis.text.y=element_text(size = 16), title = element_text(size=20)) +
  labs(title="Confronto tra densita' empirica e teorica (std)", x="rendimenti", y="Frequenza") +
  stat_function(fun=dstd,args=list(mean=fit.std$par[1],sd=fit.std$par[2],nu=fit.std$par[3]),col=2,alpha=0.7,lwd=1)
#confronto ora la funzione di ripartizione empirica con la vera funzione di ripartizione:
plot(sort(lret),(1:N/N),main="Probabilita'",col="steelblue4")
curve(pstd(x,mean=fit.std$par[1],sd=fit.std$par[2],nu=fit.std$par[3]),add=T,col=2,lwd=2)



#Punto 4.
#include.mean=F perchè la media è zero e arima(0,0,0)
spec1<-ugarchspec(variance.model=list(model="sGARCH",submodel="GARCH",garchOrder=c(1,1)),
                  mean.model=list(armaOrder=c(0,0),include.mean=F),
                  distribution.model="std")
fit1<-ugarchfit(spec1,data=lret)
show(fit1) #alpha1+beta1 =~ 1, sospetta radoce unitaria 
plot(fit1, which=9) #la distribuzone stdent-t sembra adattarsi molto bene ai nostri dati
plot(fit1, which=11) 
#dai residui standardizzati al quadrato notiamo che il primo ritardo non è significativamente 0 

spec2<-ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,2)),
                  mean.model=list(armaOrder=c(0,0),include.mean=F),
                  distribution.model="std")
fit2<-ugarchfit(spec2,data=lret) #anche in questo caso alpha1+beta1+beta2 =~ 1
#beta 2 risulta non significativo
show(fit2) 
plot(fit2, which=11)
#non si nota alcun miglioramento sostanziale rispetto al caso precedente

spec3<-ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(2,1)),
                  mean.model=list(armaOrder=c(0,0),include.mean=F),
                  distribution.model="std")
fit3<-ugarchfit(spec3,data=lret) #alpha1+alpha2+beta1 =~ 1
show(fit3)
plot(fit3, which=11)
#dall'analisi dei residui standardizzati, questo modello risulta essere migliore dei precedenti, poichè è stato in grado di ripulire la struttura di correlazione al primo lag

spec3.1<-ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(2,1)),
                  mean.model=list(armaOrder=c(0,0),include.mean=F),
                  distribution.model="std",fixed.pars=list(omega=0))
fit3.1<-ugarchfit(spec3.1,data=lret)
show(fit3.1)
plot(fit3.1, which=11)

#vista la possibile presenza di radice unitaria proviamo ad adattare un IGARCH
spec4<-ugarchspec(variance.model=list(model="iGARCH",garchOrder=c(2,1)),
                  mean.model=list(armaOrder=c(0,0),include.mean=F),distribution.model="std")
fit4<-ugarchfit(spec4,data=lret)
show(fit4)
plot(fit4, which=11)
#per quanto riguarda i residui standardizzati al quadrato, non si notano differenze significative rispetto allo standard GARCH

which.min(c(infocriteria(fit1)[1],infocriteria(fit2)[1],infocriteria(fit3)[1],infocriteria(fit3.1)[1],infocriteria(fit4)[1]))
which.min(c(infocriteria(fit1)[2],infocriteria(fit2)[2],infocriteria(fit3)[2],infocriteria(fit3.1)[2],infocriteria(fit4)[2]))
which.min(c(infocriteria(fit1)[3],infocriteria(fit2)[3],infocriteria(fit3)[3],infocriteria(fit3.1)[3],infocriteria(fit4)[3]))
which.min(c(infocriteria(fit1)[4],infocriteria(fit2)[4],infocriteria(fit3)[4],infocriteria(fit3.1)[4],infocriteria(fit4)[4]))
#numericamente il modello che minimizza tutti i criteri d'informazione è l'IGARCH(2,1) con innovazioni Student-t
#tuttavia la differenza tra questo modello e un normale GARCH(2,1) sempre con innovazioni Student-t minima
cbind(infocriteria(fit3), infocriteria(fit4))
 
persistence(fit3)
halflife(fit3) #e' molto elevato: dopo 570 giorni la varianza ritorna alla sua media
persistence(fit4)
halflife(fit4) #ovviamente l'halflife per l'igarch e' log(1) = -Inf

plot(sigma(fit3))
plot(sigma(fit4))

std.res.3=fit3@fit$residuals/sigma(fit3)
std.res.4=fit4@fit$residuals/sigma(fit4)
all.equal(std.res.3, std.res.4) #leggerissima differenza

#possiamo testare anche la presenza di effetti leverage sui redidui del modello appena stimato per vedere se è necessaria l'introduzione di un modello asimmetrico
res.std.sq<-(fit3@fit$residuals/sigma(fit3))^2
res<-fit3@fit$residuals
a=res.std.sq[2:length(res.std.sq)]
b=res[1:length(res)-1]

d=numeric(length(b))
for(i in 1:length(b))
{
  if (b[i]<0) d[i]=1 # Generazione della variabile dummy
}
summary(lm(a~1+d)) # Sign test
D=d*b
summary(lm(a~1+D)) # Size-sign test
summary(lm(a~1+d+D)) # Test congiunto
#tutti i test indicano l'assenza di effetto leverage

#Punto 5.
#previsioni della media condizionata 
n<-length(lret)
lret2<-lret[-((n-49):n)]

fit3.2<-ugarchfit(spec3,lret2)
fit4.2<-ugarchfit(spec4,lret2)
n2<-length(lret2)
prev3<-ugarchforecast(fit3.2,n.ahead=50,n.roll=0)
prev4<-ugarchforecast(fit4.2,n.ahead=50,n.roll=0)

plot(prev3,which=1)
plot(prev3,which=3)
plot(prev4,which=1)
plot(prev4,which=3)


#Punto 6.
#Calcoliamo i limiti inferiore e superiore dell'intervallo di confidenza al livello del 5%

prev<-data.frame(
  giorno<-day[(n-49):n],
  prev.garch<-prev3@forecast$seriesFor,
  prev.igarch<-prev4@forecast$seriesFor,
  lim.inf.garch<-prev3@forecast$seriesFor-2*prev3@forecast$sigmaFor,
  lim.sup.garch<-prev3@forecast$seriesFor+2*prev3@forecast$sigmaFor,
  lim.inf.igarch<-prev4@forecast$seriesFor-2*prev4@forecast$sigmaFor,
  lim.sup.igarch<-prev4@forecast$seriesFor+2*prev4@forecast$sigmaFor
)

ggplot(data = dexuseu[c(5000:n),], aes(x=day, y=lret)) + geom_line(color="steelblue", alpha=0.7) + geom_line(data=prev, aes(x=giorno, y=prev.garch), col="springgreen")+
 geom_line(data=prev, aes(x=giorno, y=lim.inf.garch), col="navyblue", lty=3)+ geom_line(data=prev, aes(x=giorno, y=lim.sup.garch), col="navyblue", lty=3) +
  geom_line(data=prev, aes(x=giorno, y=prev.igarch), col="springgreen4") + geom_line(data=prev, aes(x=giorno, y=lim.inf.igarch), col="firebrick2", lty=3)+ geom_line(data=prev, aes(x=giorno, y=lim.sup.igarch), col="firebrick2", lty=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_date(date_labels="%d %b %Y") + labs(title = "Intervalli di previsione", x = "time", y = "index")
#solitamente nei modelli garch stazionari l'intervallo di previsione dopo un certo numero di ritardi tende a stabilizzarsi. In questo caso invece, la presenza di radice pressocchè unitaria rende il comportamento del modello molto 
#smile all'igarch, il quale e' caratterizzato da una varianza esplosiva, che rende le previsioni instabili all'aumentare dell'orizzonte temporale. Dal grafico si nota infatti che entrambi gli intervalli divergono.


#Punto 7.
alpha<-c(0.01,0.05,0.10,0.90,0.95,0.99)
rollwin<-2000
n<-length(lret)
alfa.q<-qstd(1-alpha,nu=fit3@fit$coef["shape"])

#Processo GARCH: valutiamo un GARCH(2,1) senza intercetta e con distribuzione condizionata Student-t.
spec3<-ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(2,1)),
                  mean.model=list(armaOrder=c(0,0),include.mean=F),
                  distribution.model="std")
Var1.g<-Var2.g<-Var3.g<-Var4.g<-Var5.g<-Var6.g<-rep(0,n)
for(i in (rollwin+1):n) 
{
  media<-mean(lret[(i-rollwin):(i-1)])
  rend.m<-lret[i-1]-media
  mod.garch<-ugarchfit(spec=spec3,data=lret[(i-rollwin):(i-1)]-media,solver="hybrid")
  previsione<-ugarchforecast(mod.garch,n.ahead=1)
  mean<-previsione@forecast$seriesFor
  vol<-previsione@forecast$sigmaFor
  Var1.g[i]<-mean+vol*alfa.q[1]
  Var2.g[i]<-mean+vol*alfa.q[2]
  Var3.g[i]<-mean+vol*alfa.q[3]
  Var4.g[i]<-mean+vol*alfa.q[4]
  Var5.g[i]<-mean+vol*alfa.q[5]
  Var6.g[i]<-mean+vol*alfa.q[6]
  cat("obs: ",i,"\n")
}
VaR.g<-data.frame(Var1.g,Var2.g,Var3.g,Var4.g,Var5.g,Var6.g)
names(VaR.g)=c("VaR01.g","VaR05.g","VaR10.g","VaR90.g","VaR95.g","VaR99.g")
save(VaR.g,file="VaRg.RData")

##Metodo riskmetrics:
#equivalente ad un modello di IGARCH(1,1) con omega=0 e alpha=1-lambda
#lambda e' fissato a 0.94 per le serie giornaliere.

Var1.rm<-Var2.rm<-Var3.rm<-Var4.rm<-Var5.rm<-Var6.rm<-rep(0,n)
for(i in (rollwin+1):n)
{
  media<-mean(lret[(i-rollwin):(i-1)])
  rend.m<-lret[i-1]-media
  # modello i rendimenti al quadrato con il lisciamento esponenziale
  y<-lret[(i-rollwin):(i-1)] - media
  mod.hw<-HoltWinters(y^2,alpha=(1.0-0.94),gamma=FALSE,beta=FALSE,start.periods=1)
  #prevedo la volatilita' un passo in avanti
  vol.hw <-predict(mod.hw,1)
  Var1.rm[i]<-media+sqrt(vol.hw)*alfa.q[1]
  Var2.rm[i]<-media+sqrt(vol.hw)*alfa.q[2]
  Var3.rm[i]<-media+sqrt(vol.hw)*alfa.q[3]
  Var4.rm[i]<-media+sqrt(vol.hw)*alfa.q[4]
  Var5.rm[i]<-media+sqrt(vol.hw)*alfa.q[5]
  Var6.rm[i]<-media+sqrt(vol.hw)*alfa.q[6]
  cat("obs: ", i, "\n")
}
VaR.rm<-data.frame(Var1.rm,Var2.rm,Var3.rm,Var4.rm,Var5.rm,Var6.rm)
names(VaR.rm)=c("VaR01.rm","VaR05.rm","VaR10.rm","VaR90.rm","VaR95.rm","VaR99.rm")
save(VaR.rm,file="VaRrm.RData")

#Approccio bootstrap:
B<-100
Var1.b<-Var2.b<-Var3.b<-Var4.b<-Var5.b<-Var6.b<-rep(0,n)
for(i in (rollwin+1):n)
{
  mean<-mean(lret[(i-rollwin):(i-1)])
  rend.m <-lret[(i-rollwin):(i-1)]-mean
  Var1.boot<-Var2.boot<-Var3.boot<-Var4.boot<-Var5.boot<-Var6.boot<-0
  for(j in 1:B)
  {
    rend.boot<-sample(rend.m,replace=TRUE) #ricampiono i rendimenti
    Var1.boot<-Var1.boot+quantile(rend.boot,alpha[1]) #somma di 100 quantili stimati sul campione
    Var2.boot<-Var2.boot+quantile(rend.boot,alpha[2])
    Var3.boot<-Var3.boot+quantile(rend.boot,alpha[3])
    Var4.boot<-Var4.boot+quantile(rend.boot,alpha[4])
    Var5.boot<-Var5.boot+quantile(rend.boot,alpha[5])
    Var6.boot<-Var6.boot+quantile(rend.boot,alpha[6])
  }
  Var1.b[i]<-mean+Var1.boot/B #divido per B per ottenerne la media
  Var2.b[i]<-mean+Var2.boot/B
  Var3.b[i]<-mean+Var3.boot/B
  Var4.b[i]<-mean+Var4.boot/B
  Var5.b[i]<-mean+Var5.boot/B
  Var6.b[i]<-mean+Var6.boot/B
  cat("obs: ", i, "for j: ",j,"\n")
}
VaR.b<-data.frame(Var1.b,Var2.b,Var3.b,Var4.b,Var5.b,Var6.b)
names(VaR.b)=c("VaR01.b","VaR05.b","VaR10.b","VaR90.b","VaR95.b","VaR99.b")
save(VaR.b,file="VaRb.RData")

#plot a livello alpha=0.05 e alpha=0.95 per i tre metodi
ggplot()+
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=dexuseu$lret[(rollwin+1):n],color="blue"),size=0.15,alpha=0.3) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.g$VaR05.g[(rollwin+1):n], color="red"), alpha=1.4, size=0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.g$VaR95.g[(rollwin+1):n], color="yellow"), alpha = 1.4, size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR05.rm[(rollwin+1):n], color="green"), alpha = 1.4, size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR95.rm[(rollwin+1):n], color="black"), alpha = 1.4, size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.b$VaR05.b[(rollwin+1):n], color="orange"), alpha = 1.4, size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.b$VaR95.b[(rollwin+1):n], color="blue"), alpha = 1.4, size = 0.35) +
  theme(axis.text.x=element_text(angle = 0, hjust = 1, size = 16), axis.text.y=element_text(size = 16), title=element_text(size=15)) +
  scale_color_manual(name="Series", breaks=c("blue", "red", "yellow","green","black","orange","blue"),
                     values=c("blue" = "steelblue", "red" = "red", "yellow"= "yellow", "green"="forestgreen","black"="black","orange"="orange","blue"="blue"),
                     labels=c("log--returns ", "VaR05-garch", "VaR95-garch","VaR05-riskmetrics","VaR95-riskmetrics","VaR05-bootstrap","VaR95-bootstrap")) +
  theme(legend.title=element_text(color = "black", size = 14),legend.text=element_text(color="black", size=14)) +
  labs(title = "DEX US-EU (1999-01-04~2020-05-01)", subtitle= "Confronto tra diversi metodi di stima del VaR05 e del VaR95", x = " ", y = " ") +
  scale_x_date(date_labels = "%Y %b %d")
#sicuramente il modello riskmetrics è quello che coglie meglio l'andamento della serie


#backtesting:
library(GAS)
int<-(rollwin+1):n #non considero i VaR nulli
VaR.all<-data.frame(VaR.g$VaR01.g[int],VaR.rm$VaR01.rm[int],VaR.b$VaR01.b[int],VaR.g$VaR05.g[int],VaR.rm$VaR05.rm[int],VaR.b$VaR05.b[int],VaR.g$VaR10.g[int],VaR.rm$VaR10.rm[int],VaR.b$VaR10.b[int],
                    VaR.g$VaR90.g[int],VaR.rm$VaR90.rm[int],VaR.b$VaR90.b[int],VaR.g$VaR95.g[int],VaR.rm$VaR95.rm[int],VaR.b$VaR95.b[int],VaR.g$VaR99.g[int],VaR.rm$VaR99.rm[int],VaR.b$VaR99.b[int])
names(VaR.all)=c("VaR01.g","VaR01.rm","VaR01.b","VaR05.g","VaR05.rm","VaR05.b","VaR10.g","VaR10.rm","VaR10.b","VaR90.g","VaR90.rm","VaR90.b","VaR95.g","VaR95.rm","VaR95.b","VaR99.g","VaR99.rm","VaR99.b")
dim(VaR.all)
#creo la matrice che riempirò con gli output del backtestingVaR per ogni metodo di stima considerato (GARCH, Riskmetrics, Bootstrap) e per ogni livello di alpha
lret.rw<-lret[(rollwin+1):n]
out<-matrix(NA,nrow=18,ncol=7)
for(i in 1:6)
{
  for(j in 1:18)
  {
    BackTest<-BacktestVaR(lret.rw,VaR.all[,j],alpha=alpha[i],Lags=4)
    out[j,1]<-BackTest$AE
    out[j,2]<-BackTest$AD[1]
    out[j,3]<-BackTest$AD[2]
    out[j,4]<-sum(y<VaR.all[,j])
    out[j,5]<-BackTest$LRuc[2]
    out[j,6]<-BackTest$LRcc[2]
    out[j,7]<-BackTest$DQ$pvalue
  }
}
colnames(out)<-c("AE","mean AD viol.","max AD viol","viol","LRuc","LRcc","DQ")
rownames(out)<-c("VaR01.g","VaR01.rm","VaR01.b","VaR05.g","VaR05.rm","VaR05.b","VaR10.g","VaR10.rm","VaR10.b","VaR90.g","VaR90.rm","VaR90.b","VaR95.g","VaR95.rm","VaR95.b","VaR99.g","VaR99.rm","VaR99.b")

library(xtable)
xtable(out,digits=4)

##AE: Actual/Expected
#per alpha=(0.01,0.05,0.1) è migliore il metodo riskmetrics
#per alpha=(0.90,0.95,0.99) è migliore il metodo bootstrap
##AD: Absolute Deviation
##LRuc: proportion of failures. H0: alpha=alpha0
#migliori riskmetrics e bst
##LRucc: conditional coverage test. H0: indipendenza delle violazioni del VaR
#migliore rm

#il metodo di stima migliore tra quelli considerati risulta essere il metodo RiskMetrics.
#ggplot del metodo riskmetrics:
ggplot()+
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=dexuseu$lret[(rollwin+1):n],color="blue"),size=0.15,alpha=0.3) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR01.rm[(rollwin+1):n], color="red"), alpha=1.4, size=0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR05.rm[(rollwin+1):n], color = "yellow"), alpha = 1.4 , size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR10.rm[(rollwin+1):n], color = "green"), alpha = 1.4 , size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR90.rm[(rollwin+1):n], color = "black"), alpha = 1.4 , size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR95.rm[(rollwin+1):n], color = "orange"), alpha = 1.4 , size = 0.35) +
  geom_line(aes(x=dexuseu$day[(rollwin+1):n],y=VaR.rm$VaR99.rm[(rollwin+1):n], color = "blue"), alpha = 1.4 , size = 0.35) +
  theme(axis.text.x=element_text(angle = 0, hjust = 1, size = 16), axis.text.y=element_text(size = 16), title=element_text(size=17)) +
  scale_color_manual(name="Series", breaks=c("blue", "red", "yellow","green","black","orange","blue"),
                     values=c("blue" = "steelblue", "red" = "red", "yellow"= "yellow", "green"="forestgreen","black"="black","orange"="orange","blue"="blue"),
                     labels=c("log--returns ", "VaR01", "VaR05","VaR10","VaR90","VaR95","VaR99")) +
  theme(legend.title=element_text(color = "black", size = 14),legend.text=element_text(color="black", size=14)) +
  labs(title = "DEX US-EU",subtitle="dal 1999-01-04 al 2020-05-01", x = " ", y = " ") +
  scale_x_date(date_labels = "%Y %b %d")# esercizio6
