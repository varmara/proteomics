# ---
# title: "Классификация пептидов и проб"
# author: "Марина Варфоломеева"
# ---

# - Пакеты (инсталлируйте при необходимости)
# Из репозитория CRAN
# install.packages(c("dendextend", "ape", "fpc", "pvclust", "gplots", "NMF"), dependencies = TRUE)

library(RColorBrewer)


## Пример: Гребешки
expr <- read.table("data/Prot_Br_H_T.csv", header = TRUE, sep = ";", row.names = 1)
fact <- read.table("data/Prot_Br_H_T_factor.csv", header = TRUE, sep = ";", row.names = 1)

# Давайте познакомимся с данными.
# - Все ли правильно открылось?
# - Сколько экспериментальных групп? И каковы объемы выборок?
# - Есть ли пропущенные значения экспрессии?
# - Нужна ли нормализация?



# RI-plot
RIP <- function(X1, X2, main = "RI-plot", pch = 19, col = "darkgreen", lpars = list(col = "blue", lwd = 2), alpha = 0.3, xlab = "Intensity", ylab = "Ratio", ...){
  # соотношение и интенсивность
  R <- log2(rowMeans(X2) / rowMeans(X1))
  I <- log10(rowMeans(X2) * rowMeans(X1))
  # прозрачный цвет
  col_btransp <- adjustcolor(col, alpha.f = alpha)
  # график
  scatter.smooth(I, R, main = main, pch = pch, xlab = xlab, ylab = ylab, col = col_btransp, lpars = lpars, ...)
  abline(h = 0)
}

# Строим RI-plot'ы по разным факторам



## Подготовка данных к кластерному анализу
# Создаем короткие имена проб
part1 <- substr(x = fact$Oxygen, start = 0, stop = 1)
part2 <- substr(x = fact$Temperature, start = 0, stop = 2)
part3 <- rep(1:5, 6)
colnames(expr_log) <- paste(part1, part2, part3, sep = "_")

# транспонируем данные, чтобы кластеризовать пробы
texpr_log <- t(expr_log)
# матрица Евклидовых расстояний
d <- dist(x = texpr_log, method = "euclidean")


## Метод ближайшего соседа в R
hc_single <- hclust(d, method = "single")

# Дерево с помощью базовой графики
# ?plot.hclust
plot(hc_single)

# Дерево в ape
library(ape)
ph_single <- as.phylo(hc_single)
# ?plot.phylo
plot(ph_single, type = "phylogram", cex = 0.7)
axisPhylo()

# Дерево в dendextend
library(dendextend)
den_single <- as.dendrogram(hc_single)
# ?plot.dendrogram
op <- par(mar = c(4, 4, 1, 4), cex = 0.7)
plot(den_single, horiz = TRUE)

# При желании можно раскрасить лейблы
# а) Вручную
# Здесь в примере просто произвольные цвета
cols <- rainbow(30)
den_single_manual <- color_labels(dend = den_single, col = cols)
plot(den_single_manual, horiz = TRUE)

# б) При помощи функции
# Функция для превращения лейблов в цвета
get_colours <- function(dend, n_chars, palette = "Dark2"){
labs <- get_leaves_attr(dend, "label")
group <- substr(labs, start = 0, stop = n_chars)
group <- factor(group)
cols <- brewer.pal(length(levels(group)), name = palette)[group]
return(cols)
}

# Применяем функцию
cols <- get_colours(dend = den_single, n_chars = 4)
den_single_c <- color_labels(dend = den_single, col = cols)
plot(den_single_c, horiz = TRUE)


## Метод отдаленного соседа в R

## Метод невзвешенного попарного среднего в R

## Метод Варда в R



# Оценка качества кластеризации
## Кофенетическая корреляция
c_single <- cophenetic(ph_single)
cor(d, as.dist(c_single))

### Стабильность кластеров (Fang and Wang (2012))
# Нужно больше 1000 итераций в реальной жизни
library(fpc)
nsel <- nselectboot(d, B = 1000, clustermethod = hclustCBI, seed = 9646, method = "average", krange=2:11)

nsel$kopt # оптимальное число кластеров
nsel$stabk # средние значения нестабильности

# График нестабильности


## Ширина силуэта
complete3 <- cutree(hclust(d), 3)
qual3<- cluster.stats(d, complete3)
qual3$clus.avg.silwidths
mean(qual3$clus.avg.silwidths)

# Оцените ширину силуэта для 4 кластеров

## Бутстреп поддержка ветвей
# "An approximately unbiased test of phylogenetic tree selection" (Shimodaria, 2002)

library(pvclust)
# итераций должно быть 1000 и больше, здесь мало для скорости
set.seed(42)
cl_boot <- pvclust(expr_log, method.hclust = "average", nboot = 100, method.dist = "euclidean")

Дерево с величинами поддержки
plot(cl_boot)
# pvrect(cl_boot) # достоверные ветвления

seplot(cl_boot)
# seplot(cl_boot, identify = TRUE)
# print(cl_boot) # все значения

print(cl_boot, which = 16)

# Повторите с большим числом итераций


# Сопоставление деревьев: Танглграммы
set.seed(395)
untang_w <- untangle_step_rotate_2side(den_compl, den_w2, print_times = F)
# танглграмма
tanglegram(untang_w[[1]], untang_w[[2]],
           highlight_distinct_edges = FALSE,
           common_subtrees_color_lines = F,
           main = "Tanglegram",
           main_left = "Left tree",
           main_right = "Right tree",
           columns_width = c(8, 1, 8),
           margin_top = 3.2, margin_bottom = 2.5,
           margin_inner = 4, margin_outer = 0.5,
           lwd = 1.2, edge.lwd = 1.2,
           lab.cex = 1, cex_main = 1)

# Постройте танглграмму из дендрограмм, полученных методом ближайшего соседа и методом Варда.

# Тепловые карты экспрессии.
library(gplots) # для тепловых карт
## Палитры для тепловых карт
pal_green <- colorpanel(75, low = "black", mid = "darkgreen", high = "yellow")
# library(spatstat) # to convert palette to grayscale
# pal_gray <- to.grey(pal_green, weights=c(1,1,1))

dat <- as.matrix(expr_log)
heatmap.2(dat, col=pal_green, scale = "none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol = 1, keysize = 1, margins = c(8, 5))

heatmap.2(dat, col=pal_green, scale = "none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol = 1, keysize = 1, margins = c(8, 5), key.par = list(mgp = c(1.5, 0.9, 0), mar = c(3, 1, 3, 0.1), cex = 1), key.title = NA, key.xlab = NA)


# Еще один вариант тепловой карты
library(NMF)
aheatmap(dat, color = "-RdBu:256", annCol = groups)



