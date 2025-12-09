# Adds standard changes to ggplot theme

Adds standard changes to ggplot theme

## Usage

``` r
add_standardTheme(plot)
```

## Arguments

- plot:

  a ggplot object

## Value

a ggplot object

## Examples

``` r
library(ggplot2)
df <- data.frame(var1 = 1:10, var2 = sample(1:100, 10))
p <- ggplot(df, aes(x = var1, y = var2)) + geom_point()
p <- add_standardTheme(p)
```
