credentials::set_github_pat("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")


# - TRANSFORMAÇÃO DOS DADOS ----

#Abrir a base de dados no R
data <- read.csv("IPB-AB-64639.csv")
str(data)
class(data)

#Transformar dados para atribuir valores
data [data == "A_A"] <- 0
data [data == "A_B"] <- 1
data [data == "B_B"] <- 2 
data [data == "?_?"] <- NA

#Excluir a primeira coluna (coluna de identificação das amostras)
data <- data[-1]
str(data)

#Transformar base de dados de caractere para numérico
data[, c(1:ncol(data))] <- sapply(data[, c(1:ncol(data))], as.numeric)
str(data)

# - FREQUÊNCIA GENOTÍPICA ----

#Contagem de genótipos "AA", "BB" e "AB" na base de dados
results.gen <- vector("list", ncol(data)) #lista vazia

for (i in 1:ncol(data)) {
  
  pop <- length(which(data[,i] == 0)) + 
    length(which(data[,i] == 1)) + 
    length(which(data[,i] == 2))
  
  (AA <- length(which(data[,i] == 0)));AA
  (AB <- length(which(data[,i] == 1)));AB
  (BB <- length(which(data[,i] == 2)));BB
  
  (fAA <- AA/pop)
  (fAB <- AB/pop)
  (fBB <- BB/pop)
  
  results.gen[[i]] <- data.frame(Sample = colnames(data[i]),
                                 AA = AA,
                                 AB = AB,
                                 BB = BB,
                                 freqAA = fAA,
                                 freqAB = fAB,
                                 freqBB = fBB
                                 
  )
  
}

#Criar tabela com as frequências genotípicas
freq_genotipica <- do.call(rbind, results.gen)

# - FREQUÊNCIA ALÉLICA ----

#Calcular a frequência alélica por marcador: F(B) = n(B)/2N
results.allele <- vector("list", ncol(data)) #lista vazia

for (i in 1:ncol(data)){
  
  FB <- sum(data[,i], na.rm = TRUE) #NAs foram removidos
  
  Pop <- length(which(data[,i] == 0)) + 
    length(which(data[,i] == 1)) + 
    length(which(data[,i] == 2))
  
  (pB <- FB/(Pop*2))
  (qA <- 1 - pB)
  ((MAF_pB <- pB >= 0.01))
  ((MAF_qA <- qA >= 0.01))
  
  results.allele[[i]] <- data.frame(p = pB, q = qA, MAF_p = MAF_pB, MAF_q = MAF_qA) 
  
  
}

#Criar tabela com as frequências alélicas
freq_alelica <- do.call(rbind, results.allele)


##Para colocar as frquências alélicas em uma tabela, primeiramente precisa-se 
##colocar as infom de p e q em uma lista. Em seguida, as informações da lista 
##serão reunidades utilizando "rbind" de uma vez só utilizando "do.call"    

# - FILTRAGEM POR INDIVÍDUO ----
# - CALL RATE ----

##SNPs que não são “lidos” em um grande número de amostras indicam que o 
##processo de genotipagem foi falho. Assim, o call rate determina o percentual
##de SNP lidos de forma eficiente.                                              


#Calcular o call rate dos marcadores em toda a base de dados

results.cr <- vector("list", ncol(data))

for (i in 1:ncol(data)) {
  
  lidos <- length(which(data[,i] == 0)) +
    length(which(data[,i] == 1)) +
    length(which(data[,i] == 2))
  
  (cr <- ((lidos/nrow(data))*100))
  ((cr90 <- cr >= 90))
  
  results.cr[[i]] <- data.frame(Sample = colnames(data[i]),
                                Lidos = lidos,
                                call_rate = cr,
                                call_rate90 = cr90
  )
  
}

#Criar tabela com call rate

call.rate <- do.call(rbind, results.cr)

#Tabela com indivíduos com call rate >= 90

call_rate_maior_90 <- call.rate[which(call.rate$call_rate90 == TRUE),]
View(call_rate_maior_90)

# - FILTRAGEM DA BASE DE DADOS COM O CALL RATE ----

#Transpor a base de dados (transformou dataframe em uma matriz)
data.t <- t(data)
View(data.t)

#Criar vetor com apenas a coluna call.rate90 (TRUE, FALSE)
call_rate_resumido <- call.rate[,-c(1,2,3)]
call_rate_resumido <- as.vector(call_rate_resumido)

#Unir vetor call_rate_resumido (call.rate90) à matriz data.t
data_uniao <- cbind(data.t, call_rate_resumido)

#Transformar data_uniao em data.frame
data_uniao <- as.data.frame(data_uniao)

#Renomear última coluna do data.frame (coluna TRUE, FALSE)
names(data_uniao)[1000] <- "call.rate"

#Filtrar data_uniao para reter apenas marcadores com call rate > 90
data_filtrado <- data_uniao[which(data_uniao$call.rate == TRUE),]

#Excluir coluna call.rate
data_filtrado <-data_filtrado[,-1000]
str(data_filtrado)
data_filtrado_copy <- data_filtrado

# - MINIMUM ALLELE FREQUENCY ----
# - FILTRAGEM DA BASE DE DADOS COM O MAF ALELO B (p) ----

#Criar vetor com apenas a coluna MAF_p (TRUE, FALSE)
MAF_p <- freq_alelica[,-c(1,2,4)]
MAF_p <- as.vector(MAF_p)

#Unir vetor MAF_p à matriz data.t
data_uniao_MAFp <- cbind(data.t, MAF_p)

#Transformar data_uniao_MAFp em data.frame
data_uniao_MAFp <- as.data.frame(data_uniao_MAFp)
View(data_uniao_MAFp)

#Renomear última coluna do data.frame (coluna TRUE, FALSE)
names(data_uniao_MAFp)[1000] <- "MAF_p"

#Filtrar data_uniao_MAFp para reter apenas marcadores com MAF >= 0.01 para o alelo B (p)
data_filtrado_MAFp <- data_uniao_MAFp[which(data_uniao_MAFp$MAF_p == TRUE),]

#Excluir coluna MAF_p
data_filtrado_MAFp <- data_filtrado_MAFp[,-1000]
str(data_filtrado_MAFp)
data_filtrado_MAFp_copy <- data_filtrado_MAFp

# - FILTRAGEM DA BASE DE DADOS COM O MAF ALELO A (q) ----

#Criar vetor com apenas a coluna MAF_q (TRUE, FALSE)
MAF_q <- freq_alelica[,-c(1,2,3)]
MAF_q <- as.vector(MAF_q)

#Unir vetor MAF_q à matriz data.t
data_uniao_MAFq <- cbind(data.t, MAF_q)

#Transformar data_uniao em data.frame
data_uniao_MAFq <- as.data.frame(data_uniao_MAFq)

#Renomear última coluna do data.frame (coluna TRUE, FALSE)
names(data_uniao_MAFq)[1000] <- "MAF_q"

#Filtrar data_uniao_MAFq para reter apenas marcadores com MAF >= 0.01 para o alelo A (q)
data_filtrado_MAFq <- data_uniao_MAFq[which(data_uniao_MAFq$MAF_q == TRUE),]

#Excluir coluna MAF_q
data_filtrado_MAFq <- data_filtrado_MAFq[,-1000]
str(data_filtrado_MAFq)
data_filtrado_MAFq_copy <- data_filtrado_MAFq

# - CRIAR DATA FRAMEs COM AS MARCAS FILTRADAS POR CALL RATE MAFp E MAFq ----

#Criar coluna com o nome das marcas filtradas para o call rate como valores
data_filtrado_copy$row_names <- row.names(data_filtrado_copy)

#Acessar essa coluna para vizualizar as marcas
data_filtrado_copy$row_names

#Criar um vetor com apenas o nome das marcas
marcas_call_rate <- data_filtrado_copy$row_names

#Transformar esse vetor em um data frame
marcas_call_rate <- as.data.frame(marcas_call_rate)

#Inserir coluna intitulada "ordem" para posteriormente realizar o joint
marcas_call_rate$ordem <- 1:50408

#Renomear as colunas para posteriormente conseguir realizar o joint
colnames(marcas_call_rate) <- c("marcas", "ordem")
View(marcas_call_rate)

#Criar coluna com o nome das marcas filtradas para o ALELO B (p) como valores
data_filtrado_MAFp_copy$row_names <- row.names(data_filtrado_MAFp_copy)

#Acessar essa coluna para vizualizar as marcas
data_filtrado_MAFp_copy$row_names

#Criar um vetor com apenas o nome das marcas
marcas_MAFp <- data_filtrado_MAFp_copy$row_names

#Transformar esse vetor em um data frame
marcas_MAFp <- as.data.frame(marcas_MAFp)

#Inserir coluna intitulada "ordem" para posteriormente realizar o joint
marcas_MAFp$ordem1 <- 50409:101732

#Renomear as colunas para posteriormente conseguir realizar o joint
colnames(marcas_MAFp) <- c("marcas", "ordem1")
View(marcas_MAFp)

#Criar coluna com o nome das marcas filtradas para o ALELO A (q) como valores
data_filtrado_MAFq_copy$row_names <- row.names(data_filtrado_MAFq_copy)

#Acessar essa coluna para vizualizar as marcas
data_filtrado_MAFq_copy$row_names

#Criar um vetor com apenas o nome das marcas
marcas_MAFq <- data_filtrado_MAFq_copy$row_names

#Transformar esse vetor em um data frame
marcas_MAFq <- as.data.frame(marcas_MAFq)

#Inserir coluna intitulada "ordem" para posteriormente realizar o joint
marcas_MAFq$ordem2 <- 101733:149330

#Renomear as colunas para posteriormente conseguir realizar o joint
colnames(marcas_MAFq) <- c("marcas", "ordem2")
View(marcas_MAFq)

# - FILTRAGEM FINAL POR CALL RATE E MAF ----

#Filtro das marcas Call Rate e MAFp
library(dplyr)
filtro_MAFp_cr <- marcas_call_rate %>% inner_join(marcas_MAFp)
View(filtro_MAFp_cr)

#Filtro das marcas Call Rate e MAFp com MAFq
filtro_total <- filtro_MAFp_cr %>% inner_join(marcas_MAFq)
View(filtro_total)

#Criar vetor com apenas as marcas filtradas
filtro_final <- filtro_total[,-c(2,3,4)]

#Transformar vetor com as marcas filtradas em data frame
filtro_final <- as.data.frame(filtro_final)

#Renomear a coluna das marcas para posteriormente realizar o joint
colnames(filtro_final) <- c("marcas")
View(filtro_final)

#Realizar uma cópia da matriz iriginal transposta e transforma-la em data frame
data.t_copy <- as.data.frame(data.t)

#Criar uma coluna com o nome das marcas
data.t_copy$row_names <- row.names(data.t_copy)

#Visualizar a coluna das marcas do data frame inteiro
data.t_copy$row_names

#Renomear a última coluna para posteriormente realizar o joint
names(data.t_copy)[1000] <- "marcas"

#Acessar a coluna marcas do data frame inteiro
data.t_copy$marcas

#Filtrar o data frame inteiro pelo data frame que contém apenas as marcas filtradas
data_filtrado_final <- filtro_final %>% inner_join(data.t_copy)
View(data_filtrado_final)

#Fazer com que a coluna das marcas se torne rownames
rownames(data_filtrado_final) <- data_filtrado_final$marcas
View(data_filtrado_final)

#Excluir coluna que contém as marcas como valores
data_filtrado_final <- data_filtrado_final[,-1]
View(data_filtrado_final)

# - VISUALIZAÇÃO DA QUANTIDADE DE SNPS DURANTE A FILTRAGEM ----

#Criar um dataframe com todos os valores de NAs
lista_SNPs <- data.frame (Filtragem = c("Dados Brutos", "Filtragem CR", "Filtragem CR e MAFp",
                                        "Filtragem Total"),
                          Quantidade = c(64639, 
                                         50408,
                                         41841, 
                                         30306))

#Visualização da quantidade de SNPs restantes ao longo das filtragens
library(ggplot2)
ggplot(lista_SNPs, aes(y = Quantidade, x = Filtragem, fill = Quantidade)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,65000), breaks = seq(0,65000,5000)) +
  labs(title = "Quantidade de SNPs ao longo das filtragens")

# - PORCENTAGEM DE DADOS PERDIDOS POR MARCADOR ----

#Número de NAs calculado por linha (marcas) no dataframe transposto inteiro
NAs_data.t <- rowSums(is.na(data.t)) 
NAs_data.t <- as.data.frame(NAs_data.t)
head(NAs_data.t)
NAs_data.t <- colSums(NAs_data.t)
NAs_data.t <- as.numeric(NAs_data.t)
NAs_data.t
porcentagem_na_data.t <- 100
porcentagem_na_data.t

#Filtrar o data frame inteiro pelo data frame com as marcas filtradas do call rate
NAs_call_rate <- marcas_call_rate %>% inner_join(data.t_copy)

#Número de NAs calculado por linha (marcador) no dataframe filtrado pelo call rate
NAs_call_rate <- rowSums(is.na(NAs_call_rate)) 
NAs_call_rate <- as.data.frame(NAs_call_rate)
head(NAs_call_rate)
NAs_call_rate <- colSums(NAs_call_rate)
NAs_call_rate <- as.numeric(NAs_call_rate)
NAs_call_rate
porcentagem_na_call_rate <- (NAs_call_rate*100)/NAs_data.t
porcentagem_na_call_rate

#Criar dataframe com as marcas filtradas por call rate e MAFp
marcas_call_rate_MAFp <- filtro_MAFp_cr[,-c(2,3)]
marcas_call_rate_MAFp <- as.data.frame(marcas_call_rate_MAFp)
names(marcas_call_rate_MAFp)[1] <- "marcas"
View(marcas_call_rate_MAFp)

#Filtrar o data frame inteiro pelo data frame com as marcas filtradas para call rate e MAFp
NAs_MAFp <- marcas_call_rate_MAFp %>% inner_join(data.t_copy)

#Número de NAs calculado por linha no dataframe filtrado pelo call rate e MAFp
NAs_MAFp <- rowSums(is.na(NAs_MAFp)) 
NAs_MAFp <- as.data.frame(NAs_MAFp)
head(NAs_MAFp)
NAs_MAFp <- colSums(NAs_MAFp)
NAs_MAFp <- as.numeric(NAs_MAFp)
NAs_MAFp
porcentagem_na_call_rate_MAFp <- (NAs_MAFp*100)/NAs_data.t
porcentagem_na_call_rate_MAFp

#Criar dataframe com as marcas filtradas por call rate, MAFp e MAFq
marcas_call_rate_MAF_total <- filtro_total[,-c(2,3,4)]
marcas_call_rate_MAF_total <- as.data.frame(marcas_call_rate_MAF_total)
names(marcas_call_rate_MAF_total)[1] <- "marcas"
View(marcas_call_rate_MAFp)

#Filtrar o data frame inteiro pelo data frame com as marcas filtradas para call rate e MAFp
NAs_total <- marcas_call_rate_MAF_total %>% inner_join(data.t_copy)

#Número de NAs calculado por linha no dataframe filtrado pelo call rate
NAs_total <- rowSums(is.na(NAs_total)) 
NAs_total <- as.data.frame(NAs_total)
head(NAs_total)
NAs_total <- colSums(NAs_total)
NAs_total <- as.numeric(NAs_total)
NAs_total
porcentagem_total_filtrado <- (NAs_total*100)/NAs_data.t
porcentagem_total_filtrado

#Criar um dataframe com todos os valores de NAs
lista_NAs <- data.frame (Filtragem = c("Dados Brutos", "Filtragem CR", "Filtragem CR e MAFp",
                                       "Filtragem Total"),
                         Quantidade = c(100, 
                                        8.829009,
                                        7.930965, 
                                        6.178746))

View(lista_NAs)

#Visualização da quantidade de NAs ao longo das filtragens
library(ggplot2)
ggplot(lista_NAs, aes(y = Quantidade, x = Filtragem, fill = Quantidade)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,5)) +
  labs(title = "Quantidade de NAs ao longo das filtragens")

# - CONTAGEM DE DADOS PERDIDOS POR MARCADOR ----

#Número de dados perdidos (NA) por linha (marcador)
numero_na_linha <- rowSums(is.na(data_filtrado_final)) 
numero_na_linha <- as.data.frame(numero_na_linha)
View(numero_na_linha)
hist(numero_na_linha$numero_na_linha)

#Número de dados perdidos (NA) por coluna (indivíduos)
numero_na_coluna <- colSums(is.na(data_filtrado_final)) 
numero_na_coluna <- as.data.frame(numero_na_coluna)
View(numero_na_coluna)
hist(numero_na_coluna$numero_na_coluna)

# - IMPUTAÇÃO DE DADOS ----

#Visualizar base de dados filtrada
imputacao <- data_filtrado_final
imputacao <- as.matrix(imputacao)
View(imputacao)

#Realizar a soma por linha
soma <- as.data.frame(rowSums(imputacao, na.rm = TRUE)) 
View(soma)

#Renomear a coluna da soma
names(soma)[1] <- "soma"
View(soma)

#Realizar a média por linha
media <- (soma[1:nrow(soma),]/999)
media <- as.matrix(media)
View(media)

#Imputação
imputacao[1:nrow(imputacao),] <- ifelse(is.na(imputacao[1:nrow(imputacao),]), media[1:nrow(imputacao)], imputacao[1:nrow(imputacao),])
View(imputacao)

#Ajuste para aparecer até três casas decimais
imputacao <- format(round(imputacao, 3), nsmall = 3)
View(imputacao)

# - PCA SCRIPT 1 ----

#Análise de correlação
library(corrr)

#Visualização
library(ggcorrplot)

#Análise exploratória de dados multivariados
library(FactoMineR)

#Funções de visualização dos componentes principais
library(factoextra)

#Cálculo do desvio padrão das variáveis
sapply(data_filtrado_sem_na, sd)

#Normalização dos dados
data_normalized <- scale(data_filtrado_sem_na)
head(data_normalized)

#Matriz de correlação
corr_matrix <- cor(data_normalized)
#ggcorrplot(corr_matrix)

#Análise de PCA
data.pca <- princomp(corr_matrix)
summary(data.pca)                 

#Relação dos dois primeiros componentes com as colunas da base de dados
data.pca$loadings[, 1:2]

#Visulização dos componentes principais - Scree Plot: utilizado para visualizar 
#a importancia de cada componente e determinar o número de componentes a serem 
#mantidos
fviz_eig(data.pca, addlabels = TRUE) 

# Gráfico de indivíduos
fviz_pca_ind(data.pca,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             legend.title = "Representação"
)

# Gráfico das variáveis
fviz_pca_var(data.pca,
             col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     
             legend.title = "Contribuição"
)

# Gráfico das variáveis e indivíduos
fviz_pca_biplot(data.pca, repel = TRUE,
                col.var = "#2E9FDF", # cor das variáveis
                col.ind = "#696969"  # cor dos automoveis
)



# - PCA SCRIPT 2 ----

#Manipular os dados
library(tidyverse)

#PCA
library(stats)

#Funções de visualização dos componentes principais
library(factoextra)

#Cálculo do desvio padrão das variáveis
sapply(data_filtrado_sem_na, sd)

#Matriz de covariâncias
pca_cov <- prcomp(data_filtrado_sem_na)
summary(pca_cov)

#Matriz de correlação (variáveis padronizadas)
pca_corr <- prcomp(data_filtrado_sem_na, center = TRUE, scale = TRUE)
summary(pca_corr)

fviz_eig(pca_corr)

#Escores

summary(pca_corr)$x 

#Gráfico de indivíduos
fviz_pca_ind(pca_corr,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             legend.title = "Representação"
)

#Gráfico das variáveis
fviz_pca_var(pca_corr,
             col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     
             legend.title = "Contribuição"
)

#Gráfico das variáveis e indivíduos
fviz_pca_biplot(pca_corr, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969"  
)

# - HETEROZIGOSIDADE POR MARCADOR ----

#Dimensão do data frame

results.ho_marcador <- vector("list", ncol(data))

for (i in 1:ncol(data)) {
  
  (AA <- length(which(data[,i] == 0)))
  (AB <- length(which(data[,i] == 1)))
  (BB <- length(which(data[,i] == 2)))
  
  AABBAB <- sum(AA, BB, AB)
  Hobs <- ((AB*100)/AABBAB)
  
  results.ho_marcador[[i]] <- data.frame(Marcadores = colnames(data[i]),
                                         AA = AA,
                                         BB = BB,
                                         AB = AB,
                                         AABBAB = AABBAB,
                                         Hobs = Hobs
  )
  
}

ho_marcador <- do.call(rbind, results.ho_marcador)

# - HETEROZIGOSIDADE POR INDIVÍDUO ----

results.ho_individuos <- vector("list", ncol(data))

for (i in 1:nrow(data)) {
  
  (AA <- length(which(data[i,] == 0)))
  (AB <- length(which(data[i,] == 1)))
  (BB <- length(which(data[i,] == 2)))
  
  AABBAB <- sum(AA, BB, AB)
  Hobs <- ((AB*100)/AABBAB)
  
  results.ho_individuos[[i]] <- data.frame(Individuos = rownames(data[i]),
                                           AA = AA,
                                           BB = BB,
                                           AB = AB,
                                           AABBAB = AABBAB,
                                           Hobs = Hobs
  )
  
}

heterozigosidade_individuos <- which(call.rate$Call_rate95 == TRUE)
str(heterozigosidade_individuos)

heterozigosidade_individuos <- call.rate[which(call.rate$Call_rate95 == TRUE),]

# - VISUALIZAÇÃO DA HETEROZIGOSIDADE ----

library(ggplot2)

ho_temp <- as.data.frame(ho_marcador[1:6, 1:6])
str(ho_temp)

ggplot(data = ho_marcador,
       mapping = aes(Marcadores,
                     Hobs)) +
  geom_point(aes(Marcadores, Hobs)) +
  xlab("Marcadores") +
  ylab("Heterozigosidade observada") +
  ggtitle("Diversidade Genética")



ggplot(data = ho_marcador, aes(x = Hobs)) +
  geom_histogram(fill = "purple", color = "pink", bins = 10) +
  ggtitle("Diversidade Genética")



ho_marcador_count <- c(AA, BB, AB)
classe <- c("AA", "BB", "AB")
ho_marcador_count_1 <- data.frame(ho_marcador_count,classe)

ggplot(data = ho_marcador_count_1, aes(x="",y=ho_marcador_count, fill=as.factor(ho_marcador_count)))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)

matriz <- as.matrix(data)
matriz[1:10,1:10]
rownames(matriz) <- rownames(data)


