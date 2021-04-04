# title: "Методы выявления дифференциально-экспрессируемых белков"
# author: "Марина Варфоломеева"

# - [Код к этому занятию](04_differential_expression_analysis.R)
#
# - Данные о протеоме жабр гребешка _Pecten maximus_ из работы Artigaud et al. 2015
#     - [Prot_Br_H_T.csv](data/Prot_Br_H_T.csv)
#     - [Prot_Br_H_T_factor.csv](data/Prot_Br_H_T_factor.csv)

# # Тестирование статистических гипотез ====================
# # Способы выявления дифференциально экспрессируемых белков (The Good, The Bad, and The Ugly).
# # Fold change
# # t-тест
# # Модерируемый t-критерий

# ## Экспрессия белков у гребешков _Pecten maximus_  ====================

# Как работае t-тест мы будем разбирать на примере более
# полных данных об экспрессии белков у гребешков из работы
# @artigaud2015proteomic.

library(limma)
# Открываем данные экспрессии
expr <- read.table("data/Prot_Br_H_T.csv", header = TRUE, sep = ";", row.names = 1)
fact <- read.table("data/Prot_Br_H_T_factor.csv", header = TRUE, sep = ";", row.names = 1)

# Хорошо бы проверить, соответствует ли число проб в обоих файлах
dim(expr)
dim(fact)

# Есть ли пропущенные значения?
colSums(is.na(expr))
# Есть ли нули (они будут мешать при логарифмировании)
colSums(expr == 0)

# нормализуем и логарифмируем
expr_norm <- normalizeQuantiles(log2(expr + 1))

# Что известно о пробах?
str(fact)
table(fact$Oxygen, fact$Temperature)

# Отберем для этого примера только данные экспрессии и метаданные,
# относящиеся к гребешкам из 10 и 25 градусов, которые жили
# при нормальном количестве кислорода. Назовем получившиеся
# переменные `expr_subset` и `fact_subset`.
f_subset <- fact$Oxygen == "Normox" & fact$Temperature %in% c("10C", "25C")
expr_subset <- expr_norm[, f_subset]
fact_subset <- droplevels(fact[f_subset, ])

# ## t-тест в R

# Давайте сравним уровень экспрессии одного из белков

groups <- fact_subset$Temperature == "10C"
t.test(x = expr_subset[6, groups], y = expr_subset[6, !groups])

# Научимся добывать значение p-value.

t_result <- t.test(x = expr_subset[6, groups], y = expr_subset[6, !groups])



# Теперь мы готовы посчитать t-тест для каждого белка.

# 1) пишем функцию, которая считает t-test и добывает p-value
t_p_val <- function(x, f1, f2) {
  tryCatch(t.test(x = x[f1], y = x[f2])$p.value,
           error = function(e) NA)
}
# тестируем функцию
t_p_val(expr_subset[6, ], f1 = groups, f2 = !groups)

# 2) к каждой строке данных применяем наш t.test
pvals <- apply(X = expr_subset, MARGIN = 1, FUN = t_p_val,
               f1 = groups, f2 = !groups)

# В результате мы получаем вектор p-values
head(pvals)
class(pvals)

# Сколько белков, значимо меняющих экспрессию, мы нашли?
sum(pvals <= 0.05, na.rm = TRUE)
# Экспрессия каких белков различается?
ids_dif <- which(pvals <= 0.05)
rownames(expr_subset)[ids_dif]


# # Проблема множественных тестов ==========================

# ## Контроль FWER и FDR в R

p_bonf <- p.adjust(pvals, method = "bonferroni")
# было
head(pvals)
# стало
head(p_bonf)

# У скольких белков экспрессия значимо различается после поправки Бонферрони?
sum(p_bonf <= 0.05, na.rm = TRUE)

# Названия белков, экспрессия которых значимо различается после поправки Бонферрони?
names(pvals)[p_bonf <= 0.05]

# С поправкой Хольма
p_holm <- p.adjust(pvals, method = "holm")
sum(p_holm <= 0.05, na.rm = TRUE)

# После процедуры Беньямини-Хохберга (FDR)
p_bh <- p.adjust(pvals, method = "BH")
sum(p_bh <= 0.05, na.rm = TRUE)



# # Проблемы с обычным t-тестом ============================
# ## 1.t-статистика может не следовать t-распределению
# ## 2.Дисперсия экспрессии оценивается неточно на малых выборках
# ## 3.Разная дисперсия экспрессии



# # Moderated t-test (The Good) ============================


# Загружаем данные, создаем ExpressionSet
library(Biobase)

# Данные экспрессии
expr_data <- as.matrix(expr_subset)

# Данные о пробах
pheno_data <- fact_subset
pheno_metadata <- data.frame(
  labelDescription = c("Oxygen concentration", "Temperature"),
  row.names=c("Oxygen", "Temperature"))
pheno_data <- new("AnnotatedDataFrame",
                 data = pheno_data,
                 varMetadata = pheno_metadata)

# Данные о признаках (белках)
feature_data <- data.frame(Spot = rownames(expr_data))
rownames(feature_data) <- rownames(expr_data)
feature_metadata <- data.frame(
  labelDescription = c("Spot number"),
  row.names = c("Spot"))
f_data <- new("AnnotatedDataFrame",
              data = feature_data,
              varMetadata = feature_metadata)

# Данные об эксперименте
experiment_data <-
  new("MIAME",
      name = "Sebastien Artigaud et al.",
      lab = "lab",
      contact = "email@domain.com",
      title = "Proteomic responses to hypoxia at different temperatures in the great scallop (Pecten maximus).",
      abstract = "Abstract",
      other = list(notes = "partial dataset from Artigaud et al. 2015"))

# Собираем вместе
exp_set <-
  ExpressionSet(assayData = expr_data,
                phenoData = pheno_data,
                featureData = f_data,
                experimentData = experiment_data)


# Мы хотим сравнить уровень экспрессии каждого белка в
# группах, закодированных фактором `Temperature`

# Модельная матрица
X <- model.matrix(~ Temperature, pData(exp_set))
X

# линейная модель
fit <- lmFit(exp_set, design = X, method = "robust", maxit = 1000)
names(fit)

# Empirical Bayes statistics
efit <- eBayes(fit)
names(efit)

# Таблица дифференциально-экспрессируемых белков
topTable(efit, coef = 2)
numGenes <- nrow(exprs(exp_set))
full_list <- topTable(efit, number = numGenes)
# View(full_list)


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


MA_limma(efit, coef = 2, n = 3)


# Постройте и сравните графики:

# MA-plot первых 20 дифференциально экспрессируемых белков
MA_limma(efit, coef = 2, n = 20, text = F)
# MA-plot первых 20 дифференциально экспрессируемых белков, но таких, чтобы уровень экспрессии различался в 2 раза
MA_limma(efit, coef = 2, n = 20, text = F, lfc = 1)
# MA-plot первых 20 дифференциально экспрессируемых белков с уровнем экспрессии различающимся в 5 раз
MA_limma(efit, coef = 2, n = 20, text = F, lfc = log2(5))



# ## Сохраняем список всех белков в файл ===================
dir.create("results")
write.table(full_list, file = "results/pecten_diff_expression.csv", sep = "\t", quote = FALSE, col.names = NA)

# ## Добываем дифференциально-экспрессируемые белки для дальнейшей работы =======
# Первые 20 дифференциальных белков
my_list <- topTable(efit, coef = 2, n = 20)
# Фильтруем ExpressionSet
dif_exp_set <- exp_set[fData(exp_set)$Spot %in% my_list$Spot, ]

# ## Тепловая карта экспрессии дифференциальных белков ==========

library(gplots)
dat <- as.matrix(exprs(dif_exp_set))

# Короткие имена гребешков для графиков
part1 <- substr(x = pData(dif_exp_set)$Oxygen, start = 0, stop = 1)
part2 <- substr(x = pData(dif_exp_set)$Temperature, start = 0, stop = 2)
colnames(dat) <- make.unique(paste(part1, part2, sep = "_"))

# по сырым данным
pal_green <- colorpanel(75, low = "black", mid = "darkgreen", high = "yellow")
heatmap.2(dat, col = pal_green, scale = "none", key=TRUE, symkey = FALSE, density.info = "none", trace = "none", cexRow = 0.9, cexCol = 1, margins = c(4, 3), keysize = 0.8, key.par = list(mar = c(3, 0.1, 3, 0.1)))

# после дополнительной стандартизации по белкам
pal_blue_red <- colorpanel(75, low = "steelblue", mid = "black", high = "red")
heatmap.2(dat, col = pal_blue_red, scale = "row", key = TRUE, symkey = FALSE, density.info = "none", trace = "none", cexRow = 0.9, cexCol = 1, margins = c(4, 3), keysize = 0.8, key.par = list(mar = c(3, 0.1, 3, 0.1)))

# Рассмотрите, чем отличаются эти карты. Какая из них лучше
# подходит для представления результатов анализа
# дифференциальной экспрессии?

