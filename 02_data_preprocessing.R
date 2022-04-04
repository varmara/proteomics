# ---
# title: "Предварительная обработка данных"
# author: "Марина Варфоломеева"
# ---

# - Пакеты (инсталлируйте при необходимости):
# # Из репозитория CRAN
# install.packages(c("Hmisc", "RColorBrewer"))
# # С сайта Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("Biobase", "impute", "pcaMethods", "limma", "hexbin"))

# # Пример: протеом жабр гребешка _Pecten maximus_ -----

# Для работы мы будем использовать данные о протеоме жабр гребешка Pecten maximus от авторов пакета prot2D (Artigaud et al. 2013). Гребешков подвергали воздействию двух разных температур (15 и 25 градусов, по 6 гребешков в каждой группе). В этом исследовании, в общей сложности, было обнаружено 766 пятен.

# В файле “pecten.xlsx” на листе “exprs” хранятся необработанные данные интенсивностей пятен (raw volume data).
library(readxl)
pecten <- read_excel(path = "data/pecten.xlsx", sheet = "exprs")
head(pecten)

# Структура объекта pecten
str(pecten)

# Имеет смысл из номеров пятен сделать названия строк при помощи функции `rownames()` и удалить столбец с номерами пятен:
spot_names <- pecten$Spot
pecten <- as.matrix(pecten[, -1])
rownames(pecten) <- spot_names

# Данные о принадлежности гребешков к разным вариантам экспериментальной обработки, которые записаны на листе "pheno" в том же файле.
pecten.fac <- read_excel(path = "data/pecten.xlsx", sheet = "pheno")
head(pecten.fac)
str(pecten.fac)

# Сколько гребешков в каждом варианте эксперимента?
table(pecten.fac$Condition)


# Объект `pecten.fac` лучше превратить в обычный датафрейм. Переменную `Condition` лучше сделать фактором.
pecten.fac <- as.data.frame(pecten.fac)
pecten.fac$Condition <- factor(pecten.fac$Condition)


# # Импутация пропущенных значений -----

# ## Данные для демонстрации методов импутации
#
# В нашем примере пропущенных значений нет. В этом легко убедиться при помощи комбинации из нескольких функций.
colSums(is.na(pecten))
# В учебных целях "портим" данные пропущенными значениями (NA), чтобы создать пример для демонстрации работы методов импутации.

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

# Способы импутации:
# Вариант 1) Исключение переменных, в которых есть `NA`-----
f_na <- rowSums(is.na(spect)) < 1
ipect_none <- spect[f_na, ]

# Пришлось исключить очень много белков
dim(spect)
dim(ipect_none)

# Вариант 2) Замена `NA` средними значениями (mean substitution)-----

library(Hmisc) # для функции impute
ipect_mean <- t(apply(X = spect, MARGIN = 1, FUN = impute, fun = mean))

# Вариант 3) Замена средним по _k_-ближайшим соседям-----

library(impute)
# нужно, чтобы белки были в строках (как у нас)
knn_dat <- impute.knn(spect, k = 5)
# в результате импутации получился сложный объект - список
str(knn_dat)
# нам понадобится из него взять элемент data
ipect_knn <- knn_dat$data
# Теперь нет пропущенных значений
colSums(is.na(ipect_knn))


# Вариант 4) Импутация пропущенных значений при помощи байесовского анализа главных компонент-----
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

# ## Сравнение результатов импутации разными методами -----
# В данном случае, у нас есть полные исходные данные, поэтому мы можем для интереса проверить, какой из методов импутации дал наилучший результат.

# В качестве меры ошибки мы посчитаем корень из средней суммы квадратов отклонений исходных полных данных и восстановленных. Эта величина называется _root mean squared deviation_ и используется, например, для оценки качества предсказаний разных линейных моделей.
# Чем меньше значение RMSE, тем лучше.
#
# Иногда величину RMSE нормализуют --- делят либо на среднее значение, либо на диапазон значений. Полученная величина называется _normalized RMSE_(NRMSE). Нормализация позволяет сравнивать NRMSE для данных, измеренных в разных единицах.

# Функция для расчета RMSE.
RMSE <- function (act, imp, norm = FALSE){
  act <- as.matrix(act)
  imp <- as.matrix(imp)
  max_val <- max(rbind(act, imp))
  min_val <- min(rbind(act, imp))
  N <- nrow(act) * ncol(act)
  res <- sqrt(sum((act - imp)^2) / N)
  if (norm == TRUE) res <- res / (max_val - min_val)
  return(res)
}
# Пример расчета RMSE
RMSE(act = pecten, imp = ipect_mean)

# Вот значения RMSE
results <- list("Mean" = ipect_mean, "KNN" = ipect_knn, "BPCA" = ipect_bpca)
sapply(results, RMSE, act = pecten)

# И вот значения NRMSE
sapply(results, RMSE, act = pecten, norm = TRUE)

# Какой метод показал себя лучше всех?


# # Проблемы с использованием сырых данных экспрессии-----

# Сырые данные нельзя использовать для анализа по нескольким причинам.
#
# (1) Общий уровень интенсивностей пятен на разных сканах (и на разных образцах) может быть "смещен". Это может быть вызвано множеством факторов, от технических параметров до батч-эффекта.
#
# (2) Распределение интенсивностей пятен на одном геле асимметрично. Много пятен с низкой интенсивностью и несколько --- с большой.

# Создаем палитру и вектор цветов
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[pecten.fac$Condition]

# Строим боксплот, чтобы посмотреть на распределение
boxplot(pecten, outline = FALSE, col = cols, main = "Исходные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)


# Перед анализом данные сначала логарифмируют (чтобы сделать распределения интенсивностей более симметричными), затем нормализуют (чтобы сделать разные образцы более сравнимыми друг с другом).
pecten_log <- log2(pecten)

# Строим боксплот
boxplot(pecten_log, col = cols, main = "Логарифмированные\nданные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)

# После логарифмирования распределения интенсивностей стали более симметричными, но осталась на месте разница общего уровня экспрессии в разных образцах.


# # Нормализация-----
#
# Для того, чтобы выровнять форму распределений часто применяют квантильную нормализацию.
library(limma)
pecten_norm <- normalizeQuantiles(as.matrix(pecten_log))

boxplot(pecten_norm, col = cols, main = "Нормализованные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)

# После логарифмирования и нормализации распределение стало симметричным и приблизительно одинаковым во всех образцах --- с данными можно работать дальше.


# # MA-plot (RI-plot) -----
# - По оси X --- общий средний уровень (интенсивность) экспрессии во множестве образцов (= Intensity = Average)
# - По оси Y --- логарифм соотношения уровней экспрессии т.е. разница логарифмов уровней экспрессии (= Ratio = Mean)

# ## MA-plot для одной пробы против всех
limma::plotMA(pecten_norm, array = 1) # MA-plot из пакета `limma`
abline(h = c(-1, 0, 1), lty = c(2, 1, 2))

# ## MA-plot для сравнения двух групп проб
# Построим MA-plot для исходных данных, образцы из одного тритмента усредним.
X1 <- pecten_log[, 1:6]
X2 <- pecten_log[, 7:12]
X <- (rowMeans(X2) + rowMeans(X1)) / 2
Y <- rowMeans(X2) - rowMeans(X1)
scatter.smooth(x = X, y = Y, main = "Log-expression data", pch = 21, xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2))
abline(h = c(-1, 0, 1), lty = c(2, 1, 2))

# На графике исходных данных видно, (1) чем больше уровень экспрессии, тем больше разброс различий ; (2) график искривлен --- это видно по положению плотной массы точек в центре.

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
maplot(pecten_log[, 1:6], pecten_log[, 7:12], main = "Log-expression data")
maplot(pecten_norm[, 1:6], pecten_norm[, 7:12], main = "Normalized data")


# ## Боремся с оверплотингом (overplotting)-----

# У приведенных выше графиков есть неприятные свойства --- из-за того, что много точек данных 1) они накладываются друг на друга; 2) график долго рисуется, большой объем векторного файла при сохранении. Можно усовершенствовать график одним из способов.
maplot(pecten_norm[, 1:6], pecten_norm[, 7:12], main = "Normalized data\ntransparent markers", col = "darkgreen")

# 1) График с полупрозрачными точками на светлом фоне
# Генерируем полупрозрачные цвета
col_btransp <- adjustcolor("darkgreen", alpha.f = 0.6)
maplot(pecten_norm[, 1:6], pecten_norm[, 7:12], main = "Normalized data\ntransparent markers", col = col_btransp)

# 2) График с гексагональными ячейками
maplot_hex <- function(X1, X2, xbins = 30, main = "MA-plot,\nhexagonal binning", xlab = "Average log-expression", ylab = "Expression log-ratio", legend = 1, brpal = "YlGn", ...){
  library(hexbin)
  library(RColorBrewer)
  # Координаты
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  binned <- hexbin(cbind(X, Y), xbins = xbins)
  # Генерируем цвета
  ramp_ylgn <- colorRampPalette(brewer.pal(9, brpal)[-1])
  # График
  hexbin::plot(binned, colramp = ramp_ylgn, main = main, xlab = xlab, ylab = ylab, legend = legend, ...)
}

par(cex = 1)
maplot_hex(pecten_log[, 1:6], pecten_log[, 7:12], main = "Log-expression data,\nhexagonal binning", lcex = 0.5, brpal = "RdPu")

maplot_hex(pecten_norm[, 1:6], pecten_norm[, 7:12], main = "Normalized data,\nhexagonal binning")


# # Сохранение графиков в R-----

# ---Подготовка к сохранению---
## # Создаем в рабочей директории отдельную директорию для картинок, чтобы не захламлять. В данном случае, используем относительный путь.
dir.create(file.path("./figs"))

# pdf нужны размеры в дюймах
library(grid)
wid <- convertX(unit(12, "cm"), "inches")
hei <- convertY(unit(8, "cm"), "inches")
# Если не нужна точность, можете перевести размеры вручную.

# ---Сохранение графика---
# 1. Инициализируем "графический девайс".
# После выполнения этой строчки весь графический вывод
# перенаправляется на графическое устройство для pdf,
# а не на экран.
pdf("figs/f1.pdf", width = wid, height = hei, bg = "white", paper = "special", onefile = FALSE)
# 2. График. Нужно заранее проверить,
# что этот фрагмент кода работает отдельно.
op <- par(cex = 0.6)
boxplot(pecten_norm, col = cols, main = "Нормализованные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)
par(op)
# 3. Выключаем "графический девайс". После этого
# график будет записан в файл на диске,
# а графический вывод снова перенаправлен на экран.
dev.off()
# --- Конец сохранения графика ---

# можем встроить шрифты в pdf
# (если у вас установлен ghostscript)
embedFonts(file = "figs/f1.pdf", outfile = "figs/f1emb.pdf")

# png сам умеет переводить единицы длины-ширины.
png("figs/f1.png", width = 12, height = 8, units = "cm", res = 300, type = "cairo-png")
op <- par(cex = 0.6)
boxplot(pecten_norm, col = cols, main = "Нормализованные данные")
legend("topright", levels(pecten.fac$Condition), fill = pal, bty = "n", xpd = T)
par(op)
dev.off()


# # `ExpressionSet` Objects -----

# Результаты измерения интенсивности пятен на гелях обычно записываются в виде нескольких таблиц:
# - Данные об интенсивности пятен
# - Данные о пробах
# - Данные о белках
# - Данные об эксперименте в целом

# Класс `ExpressionSet` разработан специально для того, чтобы хранить данные из этих нескольких таблиц вместе.

# ## Создаем `ExpressionSet` вручную

# Давайте научимся создавать самостоятельно объекты `ExpressionSet`. Чтобы работать с этим классом объектов нам понадобятся функции из пакета `Biobase` с `Bioconductor`.

library(Biobase)

# ### `assayData` --- данные об интенсивности пятен
# В данном случае, это уже логарифмированные и нормализованные данные интенсивности пятен `pecten_norm`. Важно, чтобы это была матрица
is.matrix(pecten_norm)
assay_data <- pecten_norm

# ### `phenoData` --- данные о пробах
# Это аннотированный датафрейм (`AnnotatedDataFrame`), который состоит из двух частей: датафрейм с экспериментальными факторами, информацией о повторностях и т.п., а так же метаданные, в которых записаны расшифровки названий факторов. Подробнее см. в справке `?AnnotatedDataFrame`
# Названия строк в объекте с мета-данными должны быть такими же, как названия столбцов в матрице экспрессии.
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
#
# В данном случае, у нас нет никакой информации о белках в датасете о гребешках, хотя обычно такая информация бывает (появляется после MS/MS анализа пятен). Поэтому давайте сейчас в качестве тренировки создадим таблицу с данными о белках с единственной переменной --- номерами пятен. Имена строк в этой таблице должны совпадать с именами строк в данных об интенсивности пятен.
pecten.spots <- data.frame(Spot = rownames(pecten))
rownames(pecten.spots) <- rownames(assay_data)

feature_metadata <- data.frame(
  labelDescription = c("Spot number"),
  row.names = c("Spot"))

feature_data <- new("AnnotatedDataFrame",
              data = pecten.spots,
              varMetadata = feature_metadata)


# ### `experimentData` --- данные о самом эксперименте
# Их особенно важно добавить, если вы вдруг собираетесь делиться данными с кем-то. Данные об эксперименте для включения в `ExpressionSet` должны записаны в объект класса `MIAME`. В файле справки по этому объекту описаны названия и содержание полей (`?MIAME`). Мы заполним только некоторые для примера.
experiment_data <-
  new("MIAME",
      name = "Sebastien Artigaud et al.",
      lab = "lab",
      contact = "sebastien.artigaud@gmx.com",
      title = "Identifying differentially expressed proteins in two-dimensional electrophoresis experiments: inputs from transcriptomics statistical tools.",
      abstract="Abstract",
      other = list(notes = "dataset from prot2D package"))

# ### Собираем `ExpressionSet`
# Наконец, когда у нас есть все четыре элемента (на самом деле, достаточно минимум данных об интенсивности пятен), мы можем собрать из них `EspressionSet`.
eset <-
  ExpressionSet(assayData = assay_data,
                phenoData = pheno_data,
                featureData = feature_data,
                experimentData = experiment_data)


# ## Операции с `ExpressionSet` объектами.
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

# # Сохранение файлов с данными в R-----
#
# Можно сохранить данные на разных стадиях работы, но самое важное --- это сохранить данные после логарифмирования и нормализации, потому что они нужны для других видов анализа. Вообще все данные сохранять не нужно, ведь вы всегда можете их создать заново при помощи вашего скрипта.

write.csv(x = pecten_norm, file = "data/pecten_log2_normalized.csv")

# Чтобы потом прочитать данные
pecten_norm <- read.csv(file = "data/pecten_log2_normalized.csv", row.names = 1)

# Если вы создали ExpressionSet, то его тоже можно сохранить в файл, чтобы не пересоздавать.

save(eset, file = "data/pecten_eset.RData")

# Чтобы потом прочитать данные
load("data/pecten_eset.RData")


# # Задания для самостоятельной работы-----

# Для выполнения этих заданий вы можете использовать либо свои собственные данные, либо (уже логарифмированные) данные о протеоме сыворотки крови пациентов, страдающих разной степенью гиперплазии предстательной железы, из пакета `digeR` [@fan2009diger]:
#     - [prostate.xlsx](data/prostate.xlsx)
#     - [prostate.zip](data/prostate.zip)

# ## Задание 1-----

# Создайте искусственным образом пропущенные значения в данных. Потренируйтесь их заполнять разными способами. Пользуясь тем, что в этом датасете вам известны истинные значения экспрессии, сравните качество работы методов импутации.

# ## Задание 2-----

# Оцените распределение данных (боксплот, MA-plot). Выполните нормализацию, если это необходимо. Сохраните графики и нормализованные данные.

# ## Задание 3-----

# Создайте `ExpressionSet` и сохраните его.

