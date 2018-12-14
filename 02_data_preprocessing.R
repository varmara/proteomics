# title: "Предварительная обработка данных"
# author: "Марина Варфоломеева"

# # Пакеты (инсталлируйте при необходимости)
# # Из репозитория CRAN
# install.packages("Hmisc", "RColorBrewer")
# # С сайта Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("Biobase", "prot2D", "impute", "pcaMethods", "limma", "hexbin"))

# Пример: протеом жабр гребешка _Pecten maximus_

library(prot2D)
data(pecten)
data(pecten.fac)

dim(pecten)
dim(pecten.fac)

head(pecten)
pecten.fac

sum(is.na(pecten))
sapply(X = pecten, FUN = function(x)sum(is.na(x)))

# Импутация пропущенных значений.
# Сколько пропущенных значений?
colSums(is.na(pecten))



# Для этого примера напишем вредную функцию
# Функция, которая заполняет NA заданную пропорцию (frac) случайно расположенных ячеек в датафрейме (dfr)
spoil <- function(dfr, frac = 0.1, seed){
  # считаем число строк, столбцов
  rows <- nrow(dfr)
  cols <- ncol(dfr)
  # сколько значений нужно заменить на NA
  n_nas <- ceiling(frac * rows * cols)
  set.seed(seed)
  row_id <- sample(1:rows, size = n_nas, replace = TRUE)
  col_id <- sample(1:cols, size = n_nas, replace = TRUE)
  # заменяем на NA
  sapply(seq(n_nas), function(x){dfr[row_id[x], col_id[x]] <<- NA})
  return(dfr)
}

# "Портим" данные пропущенными значениями.
spect <- spoil(dfr = pecten, frac = 0.2, seed = 3194)


colSums(is.na(pecten))

colSums(is.na(spect))

## Исключение белков, в которых есть `NA`

n_nas <- rowSums(is.na(spect))

hist(n_nas)

f_na <- n_nas == 0

small_spect <- spect[f_na, ]

dim(small_spect)

## Замена `NA` средними значениями

library(Hmisc) # для функции impute

vec <- c(1, NA, 1:5)
mean(vec, na.rm = TRUE)
impute(vec)
impute(vec, fun = mean)

ipect_mean <- t(apply(X = spect, MARGIN = 1, FUN = impute, fun = mean))


## Замена `NA` средним по _k_-ближайшим соседям

library(impute)
# транспонируем, чтобы белки были в столбцах
trans_spect <- t(spect)
knn_dat <- impute.knn(trans_spect, k = 5)
# в результате импутации получился сложный объект - список
str(knn_dat)
# нам понадобится из него взять элемент data
ipect_knn <- t(knn_dat$data)
dim(ipect_knn)

## Импутация пропущенных значений при помощи байесовского анализа главных компонент

library(pcaMethods)
# транспонируем
trans_spect <- t(spect)
# центрируем и стандартизуем каждый столбец при помощи функции prep() из пакета pcaMethods.
scaled_spect <- prep(trans_spect, scale = "uv", center = TRUE, simple = FALSE)
# bpca
pc <- pca(scaled_spect$data, method="bpca", nPcs=2)
# восстановленные полные данные (центрированные и стандартизованные)
complete_obs <- completeObs(pc)
# возвращаем восстановленные данные в исходный масштаб
scaled_spect_complete <- prep(complete_obs, scale = scaled_spect$scale, center = scaled_spect$center, reverse = TRUE)
dim(scaled_spect_complete)
# транспонируем обратно
ipect_bpca <- t(scaled_spect_complete)
# убеждаемся, что размерность правильная
dim(ipect_bpca)


## Сравнение результатов импутации разными методами.

# В качестве меры ошибки мы посчитаем корень из средней суммы квадратов отклонений исходных полных данных и восстановленных. Эта величина называется _root mean squared deviation_

# Мы напишем функцию для рассчета RMSE

RMSE <- function(before, after, norm = FALSE){
  before <- as.matrix(before)
  after <- as.matrix(after)
  alldat <- cbind(before, after)
  N <- nrow(before) * ncol(before)
  enumerator <- sum((before - after)^2)
  res <- sqrt(enumerator / N)
  if(norm == TRUE) {
    res <- res / (min(alldat) - max(alldat))
  }
  return(res)
}



# замена средним
RMSE(pecten, ipect_mean)
RMSE(pecten, ipect_mean, norm = TRUE)

# замена средним по k-ближайшим соседям
RMSE(pecten, ipect_knn)
RMSE(pecten, ipect_knn, norm = TRUE)

# BPCA
RMSE(pecten, ipect_bpca)
RMSE(pecten, ipect_bpca, norm = TRUE)



# Нормализация и трансформация данных

# Боксплот исходных данных
# создаем палитру и вектор цветов
library(RColorBrewer)
pal <- brewer.pal(9, "Set1")
cols <- pal[pecten.fac$Condition]
# graphical parameters
# боксплот
boxplot(pecten, outline = FALSE, notch = T, col = cols, main = "Исходные данные")
legend("topright", levels(pecten.fac$Condition), fill = brewer.pal(9, "Set1"), bty = "n", xpd = T)


# Квантильная нормализация, игрушечный пример.
# матрица экспрессии
mat <- matrix(c(1, 7, 2, 10, 6, 3, 1, 4, 4, 7, 9, 2), ncol = 3)
rownames(mat) <- paste0("spot", 1:4)
colnames(mat) <- LETTERS[1:3]
mat

# матрица рангов
ranks <- apply(mat, 2, rank)
ranks

# переставляем значения экспрессии в порядке рангов

ranked_mat <- apply(mat, 2, function(x) x[order(x)])
ranked_mat

# считаем "цену" каждого ранга.

value <- rowMeans(ranked_mat)
value

# Подставляем цены рангов вместо рангов в исходную матрицу

norm_mat <- apply(ranks, 2, function(x) value[x])
norm_mat



# квантильная нормализация данных о протеоме гребешков
library(limma)
pecten_norm <- normalizeQuantiles(pecten)
boxplot(pecten_norm, outline = FALSE, boxwex = 0.7, notch = T, col = cols, main = "Нормализованные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)


# логарифмируем

pecten_log <- log2(pecten_norm)
boxplot(pecten_log, outline = FALSE, boxwex = 0.7, notch = T, col = cols, main = "Логарифмированные\nнормализованные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)


# RI-plot (MA-plot)

X1 <- pecten[, 1:6]
X2 <- pecten[, 7:12]
R <- log2(rowMeans(X2) / rowMeans(X1))
I <- log10(rowMeans(X2) * rowMeans(X1))

plot(I, R, main = "Raw data", pch = 21, xlab = "Intensity", ylab = "Ratio")
abline(h = 0)

# Что произойдет после нормализации?
X1 <- pecten_norm[, 1:6]
X2 <- pecten_norm[, 7:12]
R <- log2(rowMeans(X2) / rowMeans(X1))
I <- log10(rowMeans(X2) * rowMeans(X1))

plot(I, R, main = "Raw data", pch = 21, xlab = "Intensity", ylab = "Ratio")
abline(h = 0)



# из пакета prot2D
RIplot(pecten, n1 = 6, n2 = 6, main = "Raw data")
RIplot(pecten_norm, n1 = 6, n2 = 6, main = "Normalized data")


## Боремся с оверплотингом (overplotting)

# 1) График с полупрозрачными точками на светлом фоне
# Генерируем полупрозрачные цвета
col_btransp <- adjustcolor("darkgreen", alpha.f = 0.2)
plot(I, R, main = "Normalized data,\ntransparent markers", pch = 19, cex = 1.2, xlab = "Intensity", ylab = "Ratio", col = col_btransp)
abline(h = 0)

# 2) плотность распределения показана цветом
# Генерируем цвета
library(RColorBrewer)
ramp_ylgn <- colorRampPalette(brewer.pal(9,"YlGn")[-1])
col_density <- densCols(I, R, colramp = ramp_ylgn)
# Цвет фона, осей и пр.
op <- par(bg="black", fg="white", col.axis="white", col.lab="white", col.sub="white", col.main="white")
plot(I, R, col = col_density, pch = 19, cex = 0.5, main = 'Normalized data,\ncolour reflects density')
abline(h = 0)
par(op)

# 3) плотность распределения показана цветом, гексагональные ячейки
library(hexbin)
binned <- hexbin(cbind(I,R), xbins=30)
plot(binned, colramp = ramp_ylgn,
main='Normalized data,\nhexagonal binning', xlab = "Intensity", ylab = "Ratio", legend = 1)


# Сохранение графиков в R

# Создаем директорию для картинок, чтобы не захламлять рабочую директорию. В данном случае, используем относительный путь.
dir.create(file.path("./figs"))

# pdf нужны размеры в дюймах
library(grid)
wid <- convertX(unit(12, "cm"), "inches")
hei <- convertY(unit(8, "cm"), "inches")

pdf("figs/f1.pdf", width = wid, height = hei, bg = "white", paper = "special", onefile = FALSE)
op <- par(cex = 0.6)
plot(I, R, main = "Normalized data", pch = 19, xlab = "Intensity", ylab = "Ratio", col = col_btransp)
abline(h = 0)
par(op)
dev.off()
# можем встроить шрифты GhostScript
embedFonts(file = "figs/f1.pdf", outfile = "figs/f1emb.pdf")

# png сам умеет переводить единицы длины-ширины.
png("figs/f1.png", width = 12, height = 8, units = "cm", res = 300, type = "cairo-png")
op <- par(cex = 0.6)
plot(I, R, main = "Normalized data", pch = 19, xlab = "Intensity", ylab = "Ratio", col = col_btransp)
abline(h = 0)
par(op)
dev.off()


## Создаем `ExpressionSet` вручную

library(Biobase)

# данные об интенсивности пятен
expr_data <- as.matrix(pecten_log)

# данные о пробах
pheno_data <- pecten.fac
pheno_metadata <- data.frame(
  labelDescription = c("Experimental condition"),
  row.names=c("Condition"))
pheno_data <- new("AnnotatedDataFrame",
                 data = pheno_data,
                 varMetadata = pheno_metadata)

# данные о белках
feature_data <- data.frame(Spot = rownames(pecten_norm))
rownames(feature_data) <- rownames(expr_data)
feature_metadata <- data.frame(
  labelDescription = c("Spot number"),
  row.names = c("Spot"))
f_data <- new("AnnotatedDataFrame",
              data = feature_data,
              varMetadata = feature_metadata)

# данные о самом эксперименте
experiment_data <-
  new("MIAME",
      name="Sebastien Artigaud et al.",
      lab="lab",
      contact="email@domain.com",
      title="Identifying differentially expressed proteins in two-dimensional electrophoresis experiments: inputs from transcriptomics statistical tools.",
      abstract="Abstract",
      other=list(notes="dataset from prot2D package"))

# EspressionSet
eset <-
  ExpressionSet(assayData = expr_data,
                phenoData = pheno_data,
                featureData = f_data,
                experimentData = experiment_data)

## Операции с `ExpressionSet` объектами.
class(eset)
eset # то же самое, что print(eset)

# Извлекаем данные о пробах
pData(eset)
phenoData(eset)$Condition
phenoData(eset)

# Названия факторов
varLabels(eset)

# Информация о факторах
varMetadata(eset)
table(eset$Condition)

# Информация о белках
head(fData(eset))
fvarLabels(eset)
featureData(eset)

# Извлечение данных экспрессии
# exprs(eset) # полностью
exprs(eset)[1:5,1:3]

# Создание сабсетов
dim(eset)
sub_15 <- eset[, eset$Condition == "15C"]
dim(sub_15)


# Сохранение файлов с данными в R

write.table(pecten_log, file = "data/pecten_log2_normalized.csv", sep = "\t")

save(eset, file = "data/pecten_eset.RData")
