# Из репозитория CRAN
#install.packages(c("Hmisc", "RColorBrewer"))
# С сайта Bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(c("Biobase", "impute", "pcaMethods", "limma", "hexbin"))
#install.packages('gplots')
library(dplyr)
library(readxl)
library(impute)
library(zoo)
library(RColorBrewer)
library(limma)
library(ape)
library(dendextend)
library(pvclust)
library(Biobase)
library(gplots)

#датафрейм с качественными характеристиками эксперимента
mouse<-as.data.frame(read_excel("data/Data_Cortex_Nuclear.xls"))
#преобразовываем первый столбик и убираем из него символы после '_'
mouse$MouseID<-gsub('_15',' ',mouse$MouseID)
mouse$MouseID<-gsub('_14',' ',mouse$MouseID)
mouse$MouseID<-gsub('_13',' ',mouse$MouseID)
mouse$MouseID<-gsub('_12',' ',mouse$MouseID)
mouse$MouseID<-gsub('_11',' ',mouse$MouseID)
mouse$MouseID<-gsub('_10',' ',mouse$MouseID)
mouse$MouseID<-gsub('_9',' ',mouse$MouseID)
mouse$MouseID<-gsub('_8',' ',mouse$MouseID)
mouse$MouseID<-gsub('_7',' ',mouse$MouseID)
mouse$MouseID<-gsub('_6',' ',mouse$MouseID)
mouse$MouseID<-gsub('_5',' ',mouse$MouseID)
mouse$MouseID<-gsub('_4',' ',mouse$MouseID)
mouse$MouseID<-gsub('_3',' ',mouse$MouseID)
mouse$MouseID<-gsub('_2',' ',mouse$MouseID)
mouse$MouseID<-gsub('_1',' ',mouse$MouseID)
#преобразовываем переменные в факторы
mouse$Genotype<-as.factor(mouse$Genotype)
mouse$Treatment<-as.factor(mouse$Treatment)
mouse$Behavior<-as.factor(mouse$Behavior)
mouse$class<-as.factor(mouse$class)

#фильтруем данные в соотвествии с заданием (у нас генотип Ts65Dn)
mous_filtr<-mouse$Genotype=='Ts65Dn' 
mouse_fact<-mouse[mous_filtr,]
#убираем из дата фрейма данные об экспрессии
mouse_fac<-mouse_fact[,-c(2:78)]
# Сделаем "говорящие" этикетки.
#rownames(mouse_fac) <- make.unique(as.character(mouse_fac$Treatment))

#датафрейм с количественными характеристиками эксперимента
mouse_exp<-mouse_fact[,-c(1,79:82)]
str(mouse_exp)
mouse_exp<-as.matrix(mouse_exp)

#Замена пропущенных значений по k-ближайшим соседям
# транспонируем, чтобы белки были в столбцах (у нас и так белки в столбцах, ничего не транспонируем)
knn_data <- impute.knn(mouse_exp, k = 5)
# в результате импутации получился сложный объект - список
str(knn_data)
# нам понадобится из него взять элемент data
mouse_knn <- knn_data$data
# Теперь нет пропущенных значений
colSums(is.na(mouse_knn))

#соединяем датафреймы обратно 
mouse<-cbind.data.frame(mouse_fac,mouse_knn)
#усредняем каждые 15 экспериментов
mouse_mem<-as.data.frame(mouse%>%group_by(mouse$MouseID)%>%filter(Treatment=='Memantine')%>%summarise_all(mean))
mouse_mem<-mouse_mem[,-c(2:6)]
mouse_mem$Treatment<-'Memantine'
mouse_sal<-as.data.frame(mouse%>%group_by(mouse$MouseID)%>%filter(Treatment=='Saline')%>%summarise_all(mean))
mouse_sal<-mouse_sal[,-c(2:6)]
mouse_sal$Treatment<-'Saline'
mouse_new<-rbind(mouse_mem,mouse_sal)

#разделяем обратно
mouse_new_fac<-mouse_new[,c(1,79)]
mouse_new_pr<-mouse_new[,c(2:78)]

#транспонируем, чтобы белки были в строках, а пробы в столбцах
mouse_new_pr_t<-t(mouse_new_pr)


# Создаем палитру и вектор цветов
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[mouse_new_fac$Treatment]


# Строим боксплот, чтобы посмотреть на распределение
boxplot(mouse_new_pr_t, outline = FALSE, col = cols, main = "Исходные данные")
legend("topright", levels(mouse_new_fac$Treatment), fill = pal, bty = "n", xpd = T)
?legend
# Перед анализом данные сначала логарифмируют (чтобы сделать распределения интенсивностей более симметричными), затем нормализуют (чтобы сделать разные образцы более сравнимыми друг с другом).
mouse_log <- log2(mouse_new_pr_t)

# Строим боксплот
boxplot(mouse_log, col = cols, main = "Логарифмированные\nданные")
legend("topright", levels(mouse_new_fac$Treatment), fill = pal,  bty = "n", xpd = T)

# После логарифмирования распределения интенсивностей стали более симметричными, но осталась на месте разница общего уровня экспрессии в разных образцах.


# # Нормализация-----
#
# Для того, чтобы выровнять форму распределений часто применяют квантильную нормализацию.


mouse_norm <- normalizeQuantiles(as.matrix(mouse_log))

boxplot(mouse_norm, col = cols, main = "Нормализованные данные")
legend("topright", levels(mouse_new_fac$Treatment), fill = pal, bty = "n", xpd = T)

# # MA-plot (RI-plot) -----
# - По оси X --- общий средний уровень (интенсивность) экспрессии во множестве образцов (= Intensity = Average)
# - По оси Y --- логарифм соотношения уровней экспрессии т.е. разница логарифмов уровней экспрессии (= Ratio = Mean)


# ## MA-plot для сравнения двух групп проб
# Построим MA-plot для исходных данных, образцы из одного тритмента усредним.
X1 <- mouse_log[, 1:18]
X2 <- mouse_log[, 19:34]
X <- (rowMeans(X2) + rowMeans(X1)) / 2
Y <- rowMeans(X2) - rowMeans(X1)
scatter.smooth(x = X, y = Y, main = "Log-expression data", pch = 21, xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2))
abline(h = c(-1, 0, 1), lty = c(2, 1, 2))

# Чтобы не повторять код много раз, создадим функцию `maplot`, которая создает MA-plot для двух групп образцов.
maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  # Координаты
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  # График
  scatter.smooth(x = X, y = Y, main = main, pch = pch, xlab = xlab, ylab = ylab, lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

# После нормализации проблемы практически исчезнут.
maplot(mouse_log[, 1:18], mouse_log[, 19:34], main = "Log-expression data")
maplot(mouse_norm[,1:18 ], mouse_norm[,19:34 ], main = "Normalized data")

# Сделаем "говорящие" этикетки.
colnames(mouse_norm) <- make.unique(as.character(mouse_new_fac$Treatment))

# Чтобы строить деревья для проб, нам понадобится транспонировать исходные данные
tmouse_norm <- t(mouse_norm)
# Матрица расстояний =====================================

d <- dist(x = tmouse_norm, method = "euclidean")

hc_single <- hclust(d, method = "single")
ph_single<-as.phylo(hc_single)
plot(ph_single)
cor(d, as.dist(cophenetic(ph_single)))

# Метод отдаленного соседа
hc_compl <- hclust(d, method = "complete")
ph_compl <- as.phylo(hc_compl)
plot(ph_compl)
#кофенетическая корреляция
cor(d, as.dist(cophenetic(ph_compl)))

# Метод невзвешенного попарного среднего

hc_avg <- hclust(d, method = "average")
ph_avg <- as.phylo(hc_avg)
plot(ph_avg)
cor(d, as.dist(cophenetic(ph_avg)))


den_avg<-as.dendrogram(hc_avg)

# При желании можно раскрасить лейблы
# Вот простейшая функция, которая разбивает на группы по первым символам имени лейбла и раскрашивает в заданные цвета.

get_colours <- function(dend, n_chars, palette = "Dark2"){
  labs <- get_leaves_attr(dend, "label")
  group <- substr(labs, start = 0, stop = n_chars)
  group <- factor(group)
  cols <- brewer.pal(length(levels(group)), name = palette)[group]
  return(cols)
}

cols <- get_colours(dend = den_avg, n_chars = 2)
den_avg_c <- color_labels(dend = den_avg, col = cols)
plot(den_avg_c, horiz = TRUE)

# ## Бутстреп-поддержка ветвей =============================

# итераций должно быть nboot = 10000 и больше, здесь мало для скорости
cl_boot <- pvclust(mouse_norm, method.hclust = "average", nboot = 10000,
                   method.dist = "euclidean", iseed = 278456)
plot(cl_boot)


# # Moderated t-test (The Good) ============================


# Загружаем данные, создаем ExpressionSet


# Данные экспрессии
expr_dat <- as.matrix(mouse_norm)

# Данные о пробах (в моусе редусе фак скорее всего нужно поменять столбец с названиями строк как в экспрессии)
pheno_dat <- mouse_new_fac
pheno_metadat <- data.frame(
  labelDescription = c('mouse$MouseID','Treatment'),
  row.names=c( 'mouse$MouseID', 'Treatment'))
pheno_dat <- new("AnnotatedDataFrame",
                 data = pheno_dat,
                 varMetadata = pheno_metadat)

# Данные о признаках (белках)
feature_dat <- data.frame(Protein = rownames(expr_dat))
rownames(feature_dat) <- rownames(expr_dat)
feature_metadat <- data.frame(
  labelDescription = c("Name protein"),
  row.names = c("Protein"))
f_dat <- new("AnnotatedDataFrame",
             data = feature_dat,
             varMetadata = feature_metadat)

# Данные об эксперименте
experiment_dat <-
  new("MIAME",
      name = "Sebastien Artigaud et al.",
      lab = "lab",
      contact = "email@domain.com",
      title = "Mice Protein Expression Data Set")

# Собираем вместе
exp__set <-
  ExpressionSet(assayData = expr_dat,
                phenoData = pheno_dat,
                featureData = f_dat,
                experimentData = experiment_dat)


# Мы хотим сравнить уровень экспрессии каждого белка в
# группах, закодированных фактором `Temperature`

# Модельная матрица
X <- model.matrix(~ Treatment, pData(exp__set))
X

# линейная модель
fit <- lmFit(exp__set, design = X, method = "robust", maxit = 1000)
names(fit)

# Empirical Bayes statistics
efit <- eBayes(fit)
names(efit)

# Таблица дифференциально-экспрессируемых белков
topTable(efit, coef = 2)
numGenes <- nrow(exprs(exp__set))
full_list <- topTable(efit, number = numGenes)
View(full_list)


# ## MA-plot ===============================================

# Давайте создадим функцию, которая будет рисовать MA-plot,
# используя объект, возвращенный `eBayes()`.
#
# Аргументы:
#
# - `efit` - объект результатов `eBayes()`
# - `coef` - порядковый номер коэффициента линейной модели, для которого нужно сделать тест
# - `n` - число белков, которые нужно подписать
# - `sigif` - выделять ли дифференциальные белки?
# - `fdr` - уровень FDR коррекции
# - `lfc` - log fold change
# - `text` - подписывать ли n белков с сильнее всего различающейся экспрессией
# и т.д.

MA_limma <- function(efit, coef, n = 10, signif = TRUE, fdr = 0.05, lfc = 0, text = TRUE, cex.text = 0.8, col.text = "grey20", main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", pch = 19, pch.signif = 21, col = "darkgreen", alpha = 0.3, cex = 0.3, ...){
  # соотношение и интенсивность
  R <- efit$coefficients[, coef]
  I <- efit$Amean
  # прозрачный цвет
  col_btransp <- adjustcolor(col, alpha.f = alpha)
  # график
  plot(I, R, cex = cex, main = main, pch = pch, xlab = xlab, ylab = ylab, col = col_btransp, ...)
  abline(h = 0)
  # отмечаем дифференциально-экспрессируемые белки
  if(signif){
    sign <- p.adjust(efit$p.value[, coef], method = "BH") <= fdr
    large <- abs(efit$coefficients[, coef]) >= lfc
    points(I[sign & large], R[sign & large], cex = cex*2, col = "orange2", pch = pch.signif)
  }
  # подписываем первые n белков с сильнее всего различающейся экспрессией
  if(text){
    ord <- order(efit$lods[, coef], decreasing = TRUE)
    top_n <- ord[1:n]
    text(I[top_n], R[top_n], labels = efit$genes[top_n, ], pos = 4, cex = cex.text, col = col.text)
  }
}


MA_limma(efit, coef = 2, n = 20)

# ## Сохраняем список всех белков в файл ===================
dir.create("results")
write.table(full_list, file = "results/mouse_diff_expression.csv", sep = "\t", quote = FALSE, col.names = NA)

# ## Добываем дифференциально-экспрессируемые белки для дальнейшей работы =======
# Первые 20 дифференциальных белков
?topTable
my_list <- topTable(efit, coef = 2, n = 10)
# Фильтруем ExpressionSet
dif_exp_set <- exp__set[fData(exp__set)$Protein %in% my_list$Protein, ]



dat <- as.matrix(exprs(dif_exp_set))

# Короткие имена гребешков для графиков
part1 <- substr(x = pData(dif_exp_set)$Genotype, start = 0, stop = 1)
part2 <- substr(x = pData(dif_exp_set)$Treatment, start = 0, stop = 2)
colnames(dat) <- make.unique(paste(part1, part2, sep = "_"))

# по сырым данным
pal_green <- colorpanel(75, low = "black", mid = "darkgreen", high = "yellow")
heatmap.2(dat, col = pal_green, scale = "none", key=TRUE, symkey = FALSE, density.info = "none", trace = "none", cexRow = 0.9, cexCol = 1, margins = c(4, 3), keysize = 0.8, key.par = list(mar = c(3, 0.1, 3, 0.1)))

# после дополнительной стандартизации по белкам
pal_blue_red <- colorpanel(75, low = "steelblue", mid = "black", high = "red")
heatmap.2(dat, col = pal_blue_red, scale = "row", key = TRUE, symkey = FALSE, density.info = "none", trace = "none", cexRow = 0.9, cexCol = 1, margins = c(4, 3), keysize = 0.8, key.par = list(mar = c(3, 0.1, 3, 0.1)), hclustfun = )

