# ---
# title: "Методы выявления дифференциально-экспрессируемых белков"
# author: "Марина Варфоломеева"
# ---

# распределение логарифмов соотношений симметрично ####
sim_ratios <- function(n_max){
  #' Функция, которая возвращает соотношение двух
  #' случайных целых положительных чисел,
  #' лежащих в пределах от 0 до nmax
  a <- sample.int(n = n_max, size = 1)
  b <- sample.int(n = n_max, size = 1)
  return(a/b)
}

# Симулируем 100 000 соотношений
set.seed(932847)
simulated_ratios <- replicate(n = 100000, sim_ratios(n_max = 20000))

# Объединяем в датафрейм сырые соотношения и их логарифмы
dat_ratios <- data.frame(
  ratio = simulated_ratios,
  log_ratio = log2(simulated_ratios))

# рисуем боксплот
boxplot(dat_ratios, outline = F)
abline(h = 0, lty = 2, col = "red") # 0


## t-тест (The Bad) ####

# Открываем данные об экспрессии гребешков
expr <- read.table("data/Prot_Br_H_T.csv", header = TRUE, sep = ";", row.names = 1)
fact <- read.table("data/Prot_Br_H_T_factor.csv", header = TRUE, sep = ";", row.names = 1)

f_subset <- fact$Oxygen == "Normox" & fact$Temperature != "25C"
expr_subset <- expr[, f_subset]
fact_subset <- droplevels(fact[f_subset, ])

# нормализуем и логарифмируем, как на прошлых занятиях
library(limma)
expr_log <- log2(normalizeQuantiles(expr_subset))


# t-критерий для одного белка
groups <- fact_subset$Temperature == "10C"
tmp <- t.test(x = expr_log[1, groups], y = expr_log[1, !groups])

# как извлечь p-value?



# 1) пишем функцию, которая считает t-test и добывает p-value
t_p_val <- function(x, f1, f2) {
  tryCatch(t.test(x = x[f1], y = x[f2])$p.value,
    error = function(e) NA)
}
# тестируем функцию
t_p_val(expr_log[1, ], f1 = groups, f2 = !groups)

# 2) к каждой строке данных применяем наш t.test
pvals <- apply(X = expr_log, MARGIN = 1, FUN = t_p_val, f1 = groups, f2 = !groups)
# В результате мы получаем список p-values
head(pvals)
class(pvals)
# Его можно легко превратить в вектор
pvals <- unlist(pvals)


head(pvals)

# Сколько белков, достоверно меняющих экспрессию, мы нашли?
sum(pvals <= 0.05, na.rm = TRUE)
# Экспрессия каких белков различается?
names(pvals <= 0.05)

### Множественные сравнения ####

p_bonf <- p.adjust(pvals, method = "bonferroni")
head(p_bonf)

p_holm <- p.adjust(pvals, method = "holm")
sum(p_holm <= 0.05, na.rm = TRUE)

q_val <- p.adjust(pvals, method = "BH")
sum(q_val <= 0.05, na.rm = TRUE)

# У скольких белков экспрессия достоверно различается после поправки Бонферрони?
# Названия белков, экспрессия которых достоверно различается после поправки Бонферрони?


# Посчитайте пожалуйста самостоятельно, сколько достоверно различающихся белков будет найдено после поправки Хольма и после применения процедуры Беньямини-Хохберга.


## Moderated t-test (The Good) ####


#### Загружаем данные, создаем ExpressionSet
library(Biobase)
expr_data <- as.matrix(expr_log)
pheno_data <- fact_subset
pheno_metadata <- data.frame(
  labelDescription = c("Oxygen concentration", "Temperature"),
  row.names=c("Oxygen", "Temperature"))
pheno_data <- new("AnnotatedDataFrame",
                 data = pheno_data,
                 varMetadata = pheno_metadata)
feature_data <- data.frame(Spot = rownames(expr_log))
rownames(feature_data) <- rownames(expr_data)
feature_metadata <- data.frame(
  labelDescription = c("Spot number"),
  row.names = c("Spot"))
f_data <- new("AnnotatedDataFrame",
              data = feature_data,
              varMetadata = feature_metadata)
experiment_data <-
  new("MIAME",
      name="Sebastien Artigaud et al.",
      lab="lab",
      contact="email@domain.com",
      title="Proteomic responses to hypoxia at different temperatures in the great scallop (Pecten maximus).",
      abstract="Abstract",
      other=list(notes="partial dataset from Artigaud et al. 2014"))
exp_set <-
  ExpressionSet(assayData = expr_data,
                phenoData = pheno_data,
                featureData = f_data,
                experimentData = experiment_data)

#### Уровни фактора
groups <- pData(exp_set)$Temperature
table(groups)

#### Создаем модельную матрицу
X <- model.matrix(~ groups)

#### Подбираем линейную модель для _i_-того белка
fit <- lmFit(exp_set, design = X, method = "robust", maxit = 1000)
names(fit)

fit$coefficients[1, ]

#### Empirical Bayes statistics
efit <- eBayes(fit)
names(efit)

#### Таблица дифференциально-экспрессируемых белков
topTable(efit, coef = 2)
numGenes <- nrow(exprs(exp_set))
full_list <- topTable(efit, number = numGenes)
View(full_list)

#### RI-plot
RIP_limma <- function(efit, coef, n = 10, signif = TRUE, fdr = 0.05, lfc = 0, text = TRUE, cex.text = 0.8, col.text = "grey20", main = "RI-plot", xlab = "Intensity", ylab = "Ratio", pch = 19, pch.signif = 21, col = "darkgreen", alpha = 0.3, cex = 0.3, ...){
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

RIP_limma(efit, coef = 2, n = 20, text = F)
RIP_limma(efit, coef = 2, n = 20, text = F, lfc = 1)
RIP_limma(efit, coef = 2, n = 5)


#### Сохраняем список всех белков в файл
dir.create("results")
write.table(full_list, file = "results/pecten_diff_expression.csv", sep = "\t", quote = FALSE, col.names = NA)

#### Добываем дифференциально-экспрессируемые белки для дальнейшей работы
f_dif <- full_list$adj.P.Val <= 0.05 & abs(full_list$logFC) >= 1
# Находим имена пятен
names_dif <- full_list$Spot[f_dif]
# Находим индексы пятен в ExpressionSet
ids_dif_limma <- match(names_dif, fData(exp_set)$Spot)
# Фильтруем ExpressionSet
dif_exp_set <- exp_set[ids_dif_limma, ]
# Короткие имена гребешков
part1 <- substr(x = pData(dif_exp_set)$Oxygen, start = 0, stop = 1)
part2 <- substr(x = pData(dif_exp_set)$Temperature, start = 0, stop = 2)
short_names <- make.unique(paste(part1, part2, sep = "_"))
colnames(exprs(dif_exp_set)) <- short_names


#### Тепловая карта экспрессии дифференциальных белков
library(gplots)
dat <- as.matrix(exprs(dif_exp_set))

pal_green <- colorpanel(75, low = "black", mid = "darkgreen", high = "yellow")
heatmap.2(dat, col = pal_green, scale = "none", key=TRUE, symkey = FALSE, density.info = "none", trace = "none", cexRow = 0.9, cexCol = 1, margins = c(4, 3), keysize = 0.8, key.par = list(mar = c(3, 0.1, 3, 0.1)))

dev.off()
pal_blue_red <- colorpanel(75, low = "steelblue", mid = "black", high = "red")
heatmap.2(dat, col = pal_blue_red, scale = "row", key = TRUE, symkey = FALSE, density.info = "none", trace = "none", cexRow = 0.9, cexCol = 1, margins = c(4, 3), keysize = 0.8, key.par = list(mar = c(3, 0.1, 3, 0.1)))

