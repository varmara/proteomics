#' title: "Знакомство с R"
#' author: "Марина Варфоломеева"

#' # Основы языка R
#' ## Математические операции
2 + 3
36 / 2
7 * 4
5 ^ 2
sqrt(27)


#' ## Предупреждения и ошибки (warnings and errors)
sqrt(-27)
sqr(27)


#' ## Переменные (variables)
num_1 <- 1024 / 2
num_1
1238 * 3 -> num_2  # экзотический вариант
num_2

"это текст"
'это тоже текст'

#' ### Логические данные (`logical`)
TRUE # истина
FALSE # ложь
T
F

#' ## Встроенные в R константы

#' ### NA
NA + 2
NA * 0
NA / 0
sqrt(NA)

#' ### Inf
1 / 0

#' ### NAN
0 / 0
sqrt(-1)

#' ### NULL
23
sqrt(25)

#' #### Создание векторов из произвольных элементов
c(2, 4, 6)
c(-9.3, 0, 2.17, 21.3)


vect_num <- c(2, 4, 6, 8, 10, 12, 14, 16)
vect_num_1 <- c(1.3, 1.7, NA, 0.9, 1.6, 1.4)


c(1, 1, 5:9)
c(vect_num, vect_num_1)
c(100, vect_num)


c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE)
colours <- c("red", "orange", "yellow", "green", "blue", "violet")

LETTERS
letters
month.abb
month.name

1:10 # от одного до 10
-5:3 # от -5 до 3

rep(x = 1, times = 3) # 1 повторяется 3 раза
rep(x = "red", times = 5) # "red" повторяется 5 раз
rep(x = TRUE, times = 2) # TRUE повторяется 2 раза

rep(TRUE, 5) # TRUE повторяется 5 раз, аргументы без названий

vect_log <- c(rep(TRUE, 3), rep(FALSE, 3), rep(TRUE, 4))
vect_log

#' #### Адресация внутри векторов
vect_num # весь вектор
vect_num[1] # первый элемент
vect_num[3] # 3-й элемент

colours # весь вектор
colours[3:5] # 3-5 элемент
LETTERS[1:3]

vect_num[c(2, 4, 6)]
colours[c(1, 6)]
month.name[c(12, 1, 2)]


vect_num[2, 4, 6]
colours[1, 6]


vect_num[198]
letters[33]
month.name[13]

#' #### Векторизованные операции
vect_num + 2
vect_num * 2
vect_num * (-2)
vect_num ^2

#' ### Матрицы (matrices)

# Матрица с числовыми данными
matrix(data = 1:12, nrow = 4)
# Матрица с текстовыми данными
matrix(data = LETTERS[1:12], ncol = 6)

matrix(data = 1:6, ncol = 3)
matrix(data = 1:6, ncol = 3, byrow = TRUE)

#' #### Адресация в матрицах
mat <- matrix(data = LETTERS[1:12], ncol = 3)
mat
mat[3, 2]
mat[1, ]
mat[, 3]
mat[, -1]
mat[1:3, c(1, 3)]

#' ### Массивы (arrays)

ar <- array(data = 1:24, dim = c(2, 4, 3))
ar

#' #### Адресация в массивах
ar[1, 2, 3]
ar[, , 1]
ar[, 1, ]
ar[1, , ]
ar[1, 1:3, 1]
ar[, , -1]

#' ### Датафреймы (data frames)

x <- 1:4
y <- LETTERS[1:4]
z <- c(TRUE, TRUE, FALSE, TRUE)
dat <- data.frame(v1 = x, v2 = y, v3 = z, stringsAsFactors = FALSE)
dat

#' #### Адресация в датафреймах

dat[2, 2]

dat[, c("v1", "v3")]

dat$v1

dat$v3[3]


#' ### Списки (lists)

list(dat, mat, vect_num, colours)

lst <- list(Dfr = dat, Matr = mat, Vect1 = vect_num, Vect2 = colours)
lst

#' #### Адресация в списках

lst[1:2] # список из двух элементов
lst[1] # список из одного элемента

lst[[1]]
lst[[1]]$v2

lst$Vect1
lst$Matr[, 1]


#' ## Факторы (factors) --- особый тип данных

snail_colours <- c("red", "green", "green", "green", "yellow", "yellow", "yellow", "yellow")
snail_colours # это текстовый вектор.

f_snail_alphabet <- factor(snail_colours)
f_snail_alphabet

f_snail_ryg <- factor(snail_colours, levels = c("red", "yellow", "green"))
f_snail_ryg
f_snail_yrg <- relevel(f_snail_ryg, ref = "yellow")
f_snail_yrg

#' ## Как узнать, к какому классу структур данных относится содержимое переменной?

class(f_snail_alphabet)
class(vect_log)
class(vect_num)
class(colours)
class(mat)
class(ar)
class(dat)
class(lst)

#' ## Приведение (coercion), проверка принадлежности к классу/типу.

vect_num
as.character(vect_num)

vect_log
as.numeric(vect_log)
as.character(vect_log)
as.character(as.numeric(vect_log))

as.numeric(as.character(vect_log))

is.character(vect_num)
is.numeric(vect_num)
is.logical(vect_log)

is.vector(vect_log)
is.matrix(vect_log)
as.matrix(vect_log)

is.matrix(mat)
is.array(dat)
as.data.frame(mat)

#' ## Функции (functions)

#' Вспомним, как выглядят наши векторы
vect_num
vect_num_1

#' Длину вектора можно вычислить при помощи функции `length()`
length(vect_num)
length(vect_num_1)

#' Сумму элементов вектора при помощи функции `sum()`
sum(vect_num)
sum(vect_num_1)
?sum

#' вручную
a <- sum(vect_num_1, na.rm = TRUE) / (length(vect_num_1) - 1)
a

mean(vect_num_1)
mean(vect_num_1, na.rm = TRUE)

#' пользовательская функция
my_mean <- function(x){
  mean(x, na.rm = TRUE)
  }

my_mean <- function(x){
  res <- mean(x, na.rm = TRUE)
  return(res)
}

mean(vect_num, na.rm = TRUE)
my_mean(vect_num)

#' # Установка дополнительных пакетов c CRAN и Bioconductor
## install.packages(c("ggplot2", "gplots", "fpc", "pvclust", "Hmisc"))
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("Biobase", "prot2D", "impute", "pcaMethods", "limma"))

#' # Чтение данных из файлов

#' ## Чтение из текстовых файлов

?read.table

# https://raw.githubusercontent.com/varmara/proteomics-course/master/data/expression_3.csv

dat <- read.table(file = "data/expression_3.csv", header = TRUE, sep = ",", dec = ".")

head(dat)

sapply(dat, class)

dat1 <- within(dat, {
  X2 <- as.numeric(as.character(X2))
  X6 <- as.numeric(as.character(X6))
}
)

dat$X2[is.na(dat1$X2)]
dat$X6[is.na(dat1$X6)]

sapply(dat1, class)

#' ## Чтение из архивированных файлов

# https://raw.githubusercontent.com/varmara/proteomics-course/master/data/expression_3.zip

dat2 <- read.table(unz("./data/expression_3.zip", "expression_analysis/3.csv"), header=T, sep=",", stringsAsFactors = FALSE)
head(dat2)
sapply(dat2, class)


#' ## Чтение файлов Excel

# install.packages("readxl")
library(readxl)
dat3 <- read_excel(path = "data/expression_3.xlsx")

sapply(dat3, class)
sum(is.na(dat3$`2`))
sum(is.na(dat3$`6`))

colnames(dat3)

new_names <- paste0("X", colnames(dat3))
colnames(dat3) <- new_names

