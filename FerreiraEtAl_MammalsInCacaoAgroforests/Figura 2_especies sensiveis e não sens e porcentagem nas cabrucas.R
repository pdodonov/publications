##RRratio - 30 dias
#30 dias
#gerando Fiura 2

install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggpubr")

library("ggplot2")
library(dplyr)
library(ggplot2)
library("ggpubr")
theme_set(theme_pubclean())

ERratio_30 <- read.csv2("ERratio_30.csv")
ERratio_30
levels(ERratio_30$species)
ERratio_30$species <- ordered (ERratio_30$species, levels=c("cabtat", "sciaes","potfla","sapxan", "pectaj","pumcon","nasnas","eirbar","leochr", "daslep","didaur",
                                                            "tamtet","mazgou","cunpac","calkuh","leowie","dasnov","procan","certho"))



r <- ggplot(ERratio_30, aes(species, log10..cacao.forest., fill = ERratio_30$log10..cacao.forest. >= 0)) + geom_bar(stat = "identity")+
  labs(y = expression ('log'[10]~ 'Record ratio'), x="Species", subtitle="b)")+
  scale_color_manual(values=c("#999999","#000000"))+
  scale_fill_manual(values=c("#999999","#000000"))+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, face="italic", family="Arial"), axis.title = element_text(size = 12, family="Arial"))
r

library("RColorBrewer")
display.brewer.pal(n = 8, name = 'RdBu')
display.brewer.all()

r2 <- ggplot(ERratio_30, aes(species, log10..cacao.forest., fill = log10..cacao.forest.)) + geom_bar(stat = "identity")+
  labs(y = expression ('log'[10]~ 'Record ratio'), x="Species", subtitle="b)")+
  scale_fill_continuous(low="red",high="blue")+
  theme(panel.background = element_rect(fill = "white",colour = "white"), panel.grid.major = element_blank(),panel.grid.minor = element_line(size = 0.5, linetype = "dotted", colour = "gray"), legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, size = 11, face="italic", family="Arial"), axis.title = element_text(size = 12, family="Arial"))
r2



#colocando nomes por extenso
ERratio_30 <- read.csv2("ERratio_30nomecomp.csv")
ERratio_30
levels(ERratio_30$species)
ERratio_30$species <- ordered (ERratio_30$species, levels=c("pecaritajacu","pumaconcolor","nasuanasua","eirabarbara","leontopithecuschrysomelas", "dasyproctaleporina","didelphisaurita",
                                                            "tamanduatetradactyla","mazamagouazoubira","cuniculuspaca","callithrixkuhlii","leoparduswiedii","dasypusnovencinctus","procyoncancrivorus","cerdocyonthous"))


r.1 <- ggplot(ERratio_30, aes(species, log10..cacao.forest., fill = ERratio_30$log10..cacao.forest. >= 0)) + geom_bar(stat = "identity")+
  labs(y = expression ('log'[10]~ 'Record ratio'), x="Species")+
  scale_color_manual(values=c("#999999","#000000"))+
  scale_fill_manual(values=c("#999999","#000000"))+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, face="italic", family="Arial"), axis.title = element_text(size = 12, family="Arial"))
r.1


#fazendo o gráfico de % de sítios por região (considerando tds as agroflorestas de cacau e 100% como os 40 sítios)
#Stacked barplot with multiple groups
pr_sitios <- read.csv2("Percentagem_sitioscac.csv")
str(pr_sitios)
levels(pr_sitios$Species)
pr_sitios$Species <- ordered (pr_sitios$Species, levels=c("pectaj","pumcon","nasnas","eirbar","leochr", "daslep","didaur",
                                                          "tamtet","mazgou","cunpac","calkuh","leowie","dasnov","procan","certho"))


per <- ggplot(pr_sitios, aes(x=Species, y=Percentage, fill=Landscape)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#018571","#DFC27D"))+
  labs(y = "Percentage of sites occupied", subtitle="a)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title = element_text(size = 13, family="Arial"))+
  scale_y_continuous(limits=c(0,100))
per


#incluindo espécies que so ocorreram nas florestas tbm
pr_sitios <- read.csv2("Percentagem_sitioscaceflor.csv")
str(pr_sitios)
levels(pr_sitios$Species)
pr_sitios$Species <- ordered (pr_sitios$Species, levels=c("cabtat", "sciaes","potfla","sapxan","pectaj","pumcon","nasnas","eirbar","leochr", "daslep","didaur",
                                                          "tamtet","mazgou","cunpac","calkuh","leowie","dasnov","procan","certho"))


perf <- ggplot(pr_sitios, aes(x=Species, y=Percentage, fill=Landscape)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#018571","#DFC27D"))+
  labs(y = "Percentage of sites occupied", subtitle="a)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title = element_text(size = 12, family="Arial"))+
  scale_y_continuous(limits=c(0,100))
perf

#incluindo espécies que so ocorreram nas florestas tbm mas com porcentagens apenas para as cabrucas
pr_sitios <- read.csv2("Percentagem_sitioscaceflor0.csv")
str(pr_sitios)
levels(pr_sitios$Species)
pr_sitios$Species <- ordered (pr_sitios$Species, levels=c("cabtat", "sciaes","potfla","sapxan","pectaj","pumcon","nasnas","eirbar","leochr", "daslep","didaur",
                                                          "tamtet","mazgou","cunpac","calkuh","leowie","dasnov","procan","certho"))


perf0 <- ggplot(pr_sitios, aes(x=Species, y=Percentage, fill=Landscape)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#018571","#DFC27D"))+
  labs(y = "Percentage of sites occupied", subtitle="a)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title = element_text(size = 12, family="Arial"))+
  scale_y_continuous(limits=c(0,100))
perf0

#colocando os nomes no eixo x das especies
perf2 <- ggplot(pr_sitios, aes(x=Species, y=Percentage, fill=Landscape)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#018571","#DFC27D"))+
  labs(y = "Percentage of sites occupied",x="Species", subtitle="a)")+
  theme(legend.position="top",panel.background = element_rect(fill = "white",colour = "white"), panel.grid.major = element_blank(),panel.grid.minor = element_line(size = 0.5, linetype = "dotted", colour = "gray"), axis.text.x = element_text(angle = 60, hjust = 1, size = 11, face="italic", family="Arial"), axis.title = element_text(size = 12, family="Arial"))+
  scale_y_continuous(limits=c(0,100))
perf2



library(gridExtra)

fig <- grid.arrange(perf2,r2)
fig
ggsave(fig,filename = "Fig.2_species.tiff",dpi=600)
