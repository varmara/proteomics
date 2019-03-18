# title: "Предварительная обработка данных"
# author: "Марина Варфоломеева"

# Пакеты (инсталлируйте при необходимости): =================
# Из репозитория CRAN
install.packages(c("Hmisc", "RColorBrewer"))
# С сайта Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biobase", "impute", "pcaMethods", "limma", "hexbin"))


# # Пример: протеом жабр гребешка _Pecten maximus_ ==========

library(readxl)

# Данные экспрессии
pecten <- read_excel(path = "data/pecten.xlsx", sheet = "exprs")
head(pecten)
str(pecten)

# Имеет смысл из номеров пятен сделать названия строк при помощи функции `rownames()` и удалить столбец с номерами пятен:
spot_names <- pecten$Spot
pecten <- as.matrix(pecten[, -1])
rownames(pecten) <- spot_names

# Данные о пробах
pecten.fac <- read_excel(path = "data/pecten.xlsx", sheet = "pheno")
head(pecten.fac)
str(pecten.fac)

# Сколько гребешков в каждом варианте эксперимента?
table(pecten.fac$Condition)

# Объект `pecten.fac` лучше превратить в обычный датафрейм. Переменную `Condition` лучше сделать фактором.
pecten.fac <- data.frame(pecten.fac)
pecten.fac$Condition <- factor(pecten.fac$Condition)


# # Импутация пропущенных значений. =========================

colSums(is.na(pecten))

# ## Данные для демонстрации методов импутации

# Сколько всего чисел в pecten?
N <- prod(dim(pecten))
N

# зерно генератора случайных чисел
set.seed(392408154)

# выбираем 1000 случайных ячеек
id <- sample(1:N, size = 1000)
spect <- as.matrix(pecten)
spect[id] <- NA

# Вот что получилось:
colSums(is.na(spect))

table(rowSums(is.na(spect)))


# ## Исключение переменных, в которых есть `NA` =============

f_na <- rowSums(is.na(spect)) < 1
ipect_none <- spect[f_na, ]


# Сравним размеры получившихся датафреймов
dim(spect)
dim(ipect_none)


# ## Замена `NA` средними значениями ========================

library(Hmisc) # для функции impute

ipect_mean <- t(apply(X = spect, MARGIN = 1, FUN = impute, fun = mean))


# ## Замена средним по _k_-ближайшим соседям ================

library(impute)

# транспонируем, чтобы белки были в столбцах
trans_spect <- t(spect)
knn_dat <- impute.knn(trans_spect, k = 5)

# в результате импудации получился сложный объект - список
str(knn_dat)

# нам понадобится из него взять элемент data
ipect_knn <- t(knn_dat$data)


# ## Импутация при помощи байесовского анализа главных компонент ======

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

# ## Сравнение результатов импутации разными методами.

RMSE <- function (act, imp, norm = FALSE){
  act <- as.matrix(act)
  imp <- as.matrix(imp)
  max_val <- max(rbind(act, imp))
  min_val <- min(rbind(act, imp))
  N <- nrow(act) * ncol(act)
  res <- sqrt(sum((act - imp)^2) / N)
  if (norm == TRUE) {
    res <- res / (max_val - min_val)
    }
  return(res)
}

# На самом деле нужно повторить всю процедуру, включая генерацию `NA`, много много раз --- сделать бутстреп --- здесь мы сделаем только грубую оценку.

# Вот значения RMSE
RMSE(pecten, ipect_mean)
RMSE(pecten, ipect_knn)
RMSE(pecten, ipect_bpca)

# И вот значения NRMSE
RMSE(pecten, ipect_mean, norm = TRUE)
RMSE(pecten, ipect_knn, norm = TRUE)
RMSE(pecten, ipect_bpca, norm = TRUE)


# # Проблемы с использованием сырых данных экспрессии =========

# Создаем палитру и вектор цветов
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[pecten.fac$Condition]

# Строим боксплот, чтобы посмотреть на распределение
boxplot(pecten, outline = FALSE, col = cols, main = "Исходные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)


# # Логарифмирование ================

# Логарифмируем данные
pecten_log <- log2(pecten)

# Строим боксплот
boxplot(pecten_log, col = cols, main = "Логарифмированные\nданные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)


# # Квантильная нормализация


# Игрушечный пример

# "матрица экспрессии":
mat <- matrix(c(1, 7, 2, 10, 6, 3, 1, 4, 4, 7, 9, 2), ncol = 3)
rownames(mat) <- paste0("spot", 1:4)
colnames(mat) <- LETTERS[1:3]
mat
boxplot(mat)

# Матрица рангов
ranks <- apply(mat, 2, rank)
ranks

# Теперь нужно переставить значения в каждой из переменных в порядке, заданном их рангами. Если это сделать с переменной A (первый столбец), получится
mat[ranks[, 1], 1]

# А вот и вся ранжированная матрица
ranked_mat <- apply(mat, 2, function(x) x[order(x)])
ranked_mat

# Считаем среднее значение для каждой  из строк --- "цену" каждого ранга.
value <- rowMeans(ranked_mat)
value

# Теперь эти "цены рангов" можно подставить вместо рангов в исходную матрицу. Если это сделать с первым столбцом, то получится:
value[ranks[, 1]]

# Подставляем цены рангов вместо рангов во всю исходную матрицу:
mat_norm <- apply(ranks, 2, function(x) value[x])
mat_norm

# После нормализации форма распределения всех переменных выравнялась.
boxplot(mat_norm)


# Квантильная нормализация данных о протеоме гребешков =======

library(limma)

# Квантильная нормализация
pecten_norm <- normalizeQuantiles(as.matrix(pecten_log))

boxplot(pecten_norm, col = cols, main = "Нормализованные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)


# # MA-plot (RI-plot) ========================================

# - По оси X --- общий средний уровень (интенсивность) экспрессии во множестве образцов (= Intensity = Average).
# - По оси Y --- логарифм соотношения уровней экспрессии т.е. разница логарифмов уровней экспрессии (= Ratio = Mean)


# ## MA-plot для одной пробы против всех

plotMA(pecten_norm, array = 1) # MA-plot из пакета `limma`
abline(h = c(-1, 0, 1), lty = c(2, 1, 2))


# ## MA-plot для сравнения двух групп проб

X1 <- pecten_log[, 1:6]
X2 <- pecten_log[, 7:12]
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
  scatter.smooth(x = X, y = Y,
                 main = main, pch = pch,
                 xlab = xlab, ylab = ylab,
                 lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

maplot(pecten_log[, 1:6], pecten_log[, 7:12], main = "Log-expression data")

maplot(pecten_norm[, 1:6], pecten_norm[, 7:12], main = "Normalized data")


# ## Боремся с оверплотингом (overplotting)

# График с полупрозрачными точками на светлом фоне
# Генерируем полупрозрачные цвета
col_btransp <- adjustcolor("darkgreen", alpha.f = 0.2)
maplot(pecten_norm[, 1:6], pecten_norm[, 7:12], main = "Normalized data\ntransparent markers", col = col_btransp)

# График с гексагональными ячейками
maplot_hex <- function(X1, X2, xbins = 30, main = "MA-plot,\nhexagonal binning", xlab = "Average log-expression", ylab = "Expression log-ratio", legend = 1, ...){
  library(hexbin)
  library(RColorBrewer)
  # Координаты
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  binned <- hexbin(cbind(X, Y), xbins = xbins)
  # Генерируем цвета
  ramp_ylgn <- colorRampPalette(brewer.pal(9,"YlGn")[-1])
  # График
  hexbin::plot(binned, colramp = ramp_ylgn,
               main = main,
               xlab = xlab, ylab = ylab,
               legend = legend, ...)
}

maplot_hex(pecten_log[, 1:6], pecten_log[, 7:12], main = "Log-expression data,\nhexagonal binning")

maplot_hex(pecten_norm[, 1:6], pecten_norm[, 7:12], main = "Normalized data,\nhexagonal binning")

# # Сохранение графиков в R ================================

## # Создаем директорию для картинок, чтобы не захламлять рабочую директорию. В данном случае, используем относительный путь.
dir.create(file.path("./figs"))


# pdf нужны размеры в дюймах ----
library(grid)
wid <- convertX(unit(12, "cm"), "inches")
hei <- convertY(unit(8, "cm"), "inches")

pdf("figs/f1.pdf", width = wid, height = hei, bg = "white", paper = "special", onefile = FALSE)
op <- par(cex = 0.6)
plot(I, R, main = "Normalized data", pch = 19, xlab = "Intensity", ylab = "Ratio", col = col_btransp)
abline(h = 0)
par(op)
dev.off()
# можем встроить шрифты
embedFonts(file = "figs/f1.pdf", outfile = "figs/f1emb.pdf")


# png сам умеет переводить единицы длины-ширины. ----
png("figs/f1.png", width = 12, height = 8, units = "cm", res = 300, type = "cairo-png")
op <- par(cex = 0.6)
plot(I, R, main = "Normalized data", pch = 19, xlab = "Intensity", ylab = "Ratio", col = col_btransp)
abline(h = 0)
par(op)
dev.off()


# # `ExpressionSet` Objects =================================
#
# Результаты измерения интенсивности пятен на гелях обычно записываются в виде нескольких таблиц:
#
# - Данные об интенсивности пятен --- таблица $p \times n$, где _n_ гелей записаны в столбцах, а интенсивности _p_ белков в строках.
# - Данные о пробах --- таблица $n \times q$, в которой содержится информация о _q_ свойствах проб (об экспериментальных факторах, повторностях).
# - Данные о белках --- таблица $p \times r$, в которой описаны _r_ свойств белков (например, тривиальное название, функция).
# - Данные об эксперименте в целом --- список произвольной длины, в котором содержатся свойства эксперимента и их значения (например, информация об экспериментальном объекте, имя экспериментатора, ссылка на публикацию и т.п.).


# ## Создаем `ExpressionSet` вручную =======================
#


library(Biobase)

# ### `assayData` --- данные об интенсивности пятен
# Это матрица

is.matrix(pecten_norm)
assay_data <- pecten_norm

# ### `phenoData` --- данные о пробах
# Это аннотированный датафрейм (`AnnotatedDataFrame`)

all(rownames(pecten.fac) == colnames(pecten_norm))
# Переименовываем строки в объекте с метаданными, чтобы они назывались так же, как столбцы в матрице экспрессии.
rownames(pecten.fac) <- pecten.fac[, "Sample"]

pheno_metadata <- data.frame(
  labelDescription = c("Sample name", "Experimental condition"),
  row.names=c("Sample", "Condition"))

pheno_data <- new("AnnotatedDataFrame",
                 data = pecten.fac,
                 varMetadata = pheno_metadata)

# ### `featureData` --- данные о белках, если они есть
# Это аннотированный датафрейм (`AnnotatedDataFrame`)

pecten.spots <- data.frame(Spot = rownames(pecten))
# имена строк должны совпадать с именами строк в данных об экспрессии
rownames(pecten.spots) <- rownames(assay_data)

feature_metadata <- data.frame(
  labelDescription = c("Spot number"),
  row.names = c("Spot"))

feature_data <- new("AnnotatedDataFrame",
              data = pecten.spots,
              varMetadata = feature_metadata)


# ### `experimentData` --- данные о самом эксперименте

experiment_data <-
  new("MIAME",
      name = "Sebastien Artigaud et al.",
      lab = "lab",
      contact = "sebastien.artigaud@gmx.com",
      title = "Identifying differentially expressed proteins in two-dimensional electrophoresis experiments: inputs from transcriptomics statistical tools.",
      abstract="Abstract",
      other = list(notes = "dataset from prot2D package"))


# ### Собираем `ExpressionSet`

eset <-
  ExpressionSet(assayData = assay_data,
                phenoData = pheno_data,
                featureData = feature_data,
                experimentData = experiment_data)


# ## Операции с `ExpressionSet` объектами. ===================

class(eset)
eset # то же самое, что print(eset)

# Извлекаем данные о пробах
pData(eset)
phenoData(eset)$Condition
phenoData(eset)

# Названия факторов
varLabels(eset)

# Информация о факторах (более длинное имя, например, если есть)
varMetadata(eset)

# Можно, например, посчитать, сколько в каждой группе было образцов при помощи функции `table()`
table(eset$Condition)

# Можно извлечь информацию о белках. В данном случае ее нет - просто номера пятен.
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


# # Сохранение файлов с данными в R =========================

write.csv(x = pecten_norm, file = "data/pecten_log2_normalized.csv")

# # Чтобы потом прочитать данные
# pecten_norm <- read.csv(file = "data/pecten_log2_normalized.csv", row.names = 1)


save(eset, file = "data/pecten_eset.RData")

# # Чтобы потом прочитать данные
# load("data/pecten_eset.RData")


# # Задания для самостоятельной работы ======================
#
# Для выполнения этих заданий вы можете использовать либо свои собственные данные, либо (уже логарифмированные) данные о протеоме сыворотки крови пациентов, страдающих разной степенью гиперплазии предстательной железы, из пакета `digeR` [@fan2009diger]:
#   - prostate.xlsx или prostate.zip


# ## Задание 1 ----------------------------------------------
#
# Создайте искуственным образом пропущенные значения в данных. Потренируйтесь их заполнять разными способами. Пользуясь тем, что в этом датасете вам известны истинные значения экспрессии, сравните качество работы методов импутации.



# ## Задание 2 ----------------------------------------------
#
# Оцените распределение данных (боксплот, MA-plot). Выполните нормализацию, если это необходимо. Сохраните графики и нормализованные данные.



# ## Задание 3 ----------------------------------------------
#
# Создайте `ExpressionSet` и сохраните его.


