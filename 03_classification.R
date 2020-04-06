# title: "Классификация и ординация"
# author: "Марина Варфоломеева"


# - Пакеты (инсталлируйте при необходимости) -----------------------------
## # Из репозитория CRAN
## install.packages(c("dendextend", "ape", "vegan", "pvclust", "gplots", "NMF"), dependencies = TRUE)


# ## Кластерный анализ в R: гребешки

# Вспомним, на чем мы остановились в прошлый раз.
library(readxl)
library(limma)

# Данные об экспрессии
pecten <- read_excel(path = "data/pecten.xlsx", sheet = "exprs")
spot_names <- pecten$Spot
pecten <- as.matrix(pecten[, -1])
rownames(pecten) <- spot_names

# Данные о пробах
pecten.fac <- read_excel(path = "data/pecten.xlsx", sheet = "pheno")
pecten.fac <- data.frame(pecten.fac)
pecten.fac$Condition <- factor(pecten.fac$Condition)

# Логарифмируем данные
pecten_log <- log2(pecten)
# Квантильная нормализация
pecten_norm <- normalizeQuantiles(as.matrix(pecten_log))


# Названия проб в этом файле --- длинные непонятные аббревиатуры.
colnames(pecten_norm)

colnames(pecten_norm) <- make.unique(as.character(pecten.fac$Condition))

# Чтобы строить деревья для проб, нам понадобится транспонировать исходные данные
tpecten_norm <- t(pecten_norm)


# Матрица расстояний.
d <- dist(x = tpecten_norm, method = "euclidean")

# Метод ближайшего соседа
hc_single <- hclust(d, method = "single")



# Деревья можно визуализировать при помощи базовой графики
# ?plot.hclust
plot(hc_single)

# Визуализируем средствами пакета `ape` `r citep(citation("ape"))`.
library(ape)
ph_single <- as.phylo(hc_single)
# ?plot.phylo
plot(ph_single, type = "phylogram", cex = 0.7)
axisPhylo()

# Визуализируем средствами `dendextend` `r citep(citation("dendextend"))`.
library(dendextend)
den_single <- as.dendrogram(hc_single)
# ?plot.dendrogram
op <- par(mar = c(4, 4, 1, 4), cex = 0.7)
plot(den_single, horiz = TRUE)


# При желании можно раскрасить лейблы.
library(RColorBrewer)
get_colours <- function(dend, n_chars, palette = "Dark2"){
labs <- get_leaves_attr(dend, "label")
group <- substr(labs, start = 0, stop = n_chars)
group <- factor(group)
cols <- brewer.pal(length(levels(group)), name = palette)[group]
return(cols)
}

cols <- get_colours(dend = den_single, n_chars = 2)
den_single_c <- color_labels(dend = den_single, col = cols)
plot(den_single_c, horiz = TRUE)


# ### Задание 1 --------------------------------------------
#
# Постройте дендрограммы, описывающие сходство проб, при помощи методов отдаленного соседа и среднегруппового расстояния.
hc_compl <- hclust(d, method = "complete")
ph_compl <- as.phylo(hc_compl)
plot(ph_compl)

hc_avg <- hclust(d, method = "average")
ph_avg <- as.phylo(hc_avg)
plot(ph_avg)

par(mfrow = c(1, 3), cex = 1.1)
plot(ph_single, main = "single")
plot(ph_compl, main = "complete")
plot(ph_avg, main = "average")
par(mfrow = c(1, 1), cex = 0.9)

# ## Кофенетическая корреляция =============================


# Кофенетическое расстояние
cophenetic(ph_single)


# Кофенетическая корреляция
cor(d, as.dist(cophenetic(ph_single)))


# ### Задание 2 ---------------------------------------------
#
# Оцените при помощи кофенетической корреляции качество кластеризаций, полученных разными методами. Какой метод дает лучший результат?
cor(d, as.dist(cophenetic(ph_compl)))
cor(d, as.dist(cophenetic(ph_avg)))



# ## Бутстреп-поддержка ветвей =============================


# "An approximately unbiased test of phylogenetic tree selection" (Shimodaria, 2002)

library(pvclust)
# итераций должно быть 10000 и больше, здесь мало для скорости
cl_boot <- pvclust(pecten_norm, method.hclust = "average", nboot = 100,
                   method.dist = "euclidean", iseed = 278456)

plot(cl_boot)
# pvrect(cl_boot) # достоверные ветвления

# Но точно ли мы оценили AU p-values?
seplot(cl_boot)
seplot(cl_boot, identify = TRUE)
print(cl_boot)
# 0.073 для кластера 8

# ### Задание 3 --------------------------------------------
#
# Повторите бутстреп с 1000 итераций. Чему теперь будет равна стандартная ошибка AU p-value для 8 кластера. Используйте тот же сид, что в прошлом примере.
cl_boot1 <- pvclust(pecten_norm, method.hclust = "average", nboot = 1000, method.dist = "euclidean", iseed = 278456)

plot(cl_boot1)
# pvrect(cl_boot1) # достоверные ветвления

# Но точно ли мы оценили AU p-values?
seplot(cl_boot1)
# seplot(cl_boot1, identify = TRUE)
print(cl_boot1)
# 0.016



# # Танглграмма  ===========================================

set.seed(395)
untang_w <- untangle_step_rotate_2side(den_single, den_avg, print_times = F)
# танглграмма
tanglegram(untang_w[[1]], untang_w[[2]],
           highlight_distinct_edges = FALSE,
           common_subtrees_color_lines = F,
           main = "Tanglegram",
           main_left = "Single linkage",
           main_right = "UPGMA",
           columns_width = c(8, 1, 8),
           margin_top = 3.2, margin_bottom = 2.5,
           margin_inner = 4, margin_outer = 0.5,
           lwd = 1.2, edge.lwd = 1.2,
           lab.cex = 1, cex_main = 1)

# # Тепловая карта

library(gplots) # для тепловых карт

# Палитры для тепловых карт
pal_green <- colorpanel(75, low = "black", mid = "darkgreen", high = "yellow")
# library(spatstat) # to convert palette to grayscale
# pal_gray <- to.grey(pal_green, weights=c(1,1,1))

dat <- as.matrix(pecten_norm)
heatmap.2(dat, col = pal_green, scale = "none",
          key = TRUE, symkey = FALSE,
          density.info = "none", trace = "none",
          cexRow = 1, cexCol = 1,
          keysize = 1, margins = c(8, 5))

# Настройка внешнего вида
heatmap.2(dat, col = pal_green, scale = "none",
          key = TRUE, symkey = FALSE,
          density.info = "none", trace = "none",
          cexRow = 1, cexCol = 1, keysize = 1,
          margins = c(8, 5),
          key.par = list(mgp = c(1.5, 0.9, 0),
                         mar = c(3, 1, 3, 0.1), cex = 1),
          key.title = NA, key.xlab = NA)

# Еще один вариант
library(NMF)
aheatmap(dat, color = "-RdBu:256", scale = "none",
         annCol = pecten.fac$Group,
         hclustfun = "average")

# # Ординация

# ## nMDS ординация в R: гребешки ==========================

library(vegan)
pecten_ord <- metaMDS(tpecten_norm,
                      distance = "euclidean",
                      autotransform = FALSE)

# Простейший график ординации проб:
ordiplot(pecten_ord, type = "t", display = "sites")

# Хорошая ли получилась ординация, можно узнать по величине стресса.
# Добудьте ее из объекта pecten_ord



# Раскрасим график ординации.

# Палитры
pal_col <- c("steelblue", "orangered")
pal_sh <- c(17, 19)

# Украшенный график nMDS ординации
ordiplot(pecten_ord, type = "n", display = "sites")
points(pecten_ord,
       col = pal_col[pecten.fac$Condition],
       pch = pal_sh[pecten.fac$Condition])
legend("topleft",
       levels(pecten.fac$Condition),
       col = pal_col,
       pch = pal_sh,
       bty = "n",
       xpd = T)



# График nMDS ординации, где обведено облако проб одной категории

ordiplot(pecten_ord, type = "n", display = "sites")
points(pecten_ord,
       col = pal_col[pecten.fac$Condition],
       pch = pal_sh[pecten.fac$Condition])
ordihull(pecten_ord, groups = pecten.fac$Condition, col = pal_col, label = TRUE)


# График nMDS ординации с наложенной дендрограммой

ordiplot(pecten_ord, type = "n", display = "sites")
points(pecten_ord,
       col = pal_col[pecten.fac$Condition],
       pch = pal_sh[pecten.fac$Condition])
ordicluster(pecten_ord, cluster = hc_avg)
legend("topleft",
       levels(pecten.fac$Condition),
       col = pal_col,
       pch = pal_sh,
       bty = "n",
       xpd = T)


# # Задание для самостоятельной работы ---------------------
#
# Для выполнения этого задания вы можете использовать либо свои собственные данные, либо (уже логарифмированные) данные о протеоме сыворотки крови пациентов, страдающих разной степенью гиперплазии предстательной железы, из пакета `digeR` [@fan2009diger]:
#     - [prostate.xlsx](data/prostate.xlsx)
#     - [prostate.zip](data/prostate.zip)
#
# В качестве исходных данных используйте матрицу евклидовых расстояний между пробами.
#
# - Постройте дендрограмму. Используйте алгоритм кластеризации, который лучше всего отражает матрицу исходных расстояний на дендрограмме)
#
# - Постройте танглграмму из двух дендрограмм, полученных методом ближайшего соседа и методом невзвешенного попарного среднего.
#
# - Постройте тепловую карту.
#
# - Постройте ординацию методом nMDS.
#
# На ваш взгляд, для каких целей лучше всего подходит каждый из использованных методов визуализации?




