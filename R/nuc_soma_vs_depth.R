library(ggplot2)

nuc.soma <- read.csv("C:/Users/trygveb/Dropbox/AIBS/Transcriptomics/Manuscripts/WholeCell_vs_Nuc/Figures/Figure 5/nuc_soma_exp/Mouse_V1_soma_nuc_areas.csv")

nuc.ratio.df <- data.frame()
for (i in seq(1, nrow(nuc.soma), 2)) {
  x.mean <- mean(nuc.soma$CenterX[c(i, i + 1)])
  y.mean <- mean(nuc.soma$CenterY[c(i, i + 1)])
  # Check that nuc/soma outlines are paired (centers aligned)
  x.cv <- sd(nuc.soma$CenterX[c(i, i + 1)]) / mean(nuc.soma$CenterX[c(i, i + 1)])
  nuc.area <- nuc.soma$Area[i + 1]
  soma.area <- nuc.soma$Area[i]
  nuc.area.ratio <- nuc.area / soma.area
  nuc.vol.ratio <- nuc.area.ratio^(3/2)
  nuc.ratio.df <- rbind(nuc.ratio.df, cbind(x.mean, y.mean, x.cv, nuc.area, soma.area,
                                            nuc.area.ratio, nuc.vol.ratio))
}

# Rotate coords
alpha <- -68  # Measured in photoshop
rotm <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol=2)
xy <- cbind(nuc.ratio.df$x.mean, -nuc.ratio.df$y.mean)
xy.rot <- t(rotm %*% (t(xy)-c(xy[1,1],xy[1,2]))+c(xy[1,1],xy[1,2]))
xy.rot <- sweep(xy.rot, 2, apply(xy.rot, 2, min), "-")

nuc.ratio.df <- cbind(nuc.ratio.df, xy.rot)
nuc.ratio.df$depth_from_pia <- max(nuc.ratio.df$`2`) - nuc.ratio.df$`2`
nuc.ratio.df$soma.area.bin <- cut(nuc.ratio.df$soma.area, breaks = c(0, 80, 100, 150, 250))


par(mfrow = c(2, 3))
plot(xy, asp = 1)
plot(xy.rot, asp = 1)

plot(nuc.ratio.df$depth_from_pia, nuc.ratio.df$nuc.area.ratio, 
     cex = 3 * nuc.ratio.df$soma.area,
     main = "Nuc area ratio")
lm1 <- lm(nuc.area.ratio ~ depth_from_pia, nuc.ratio.df)
summary(lm1)
abline(lm1)

soma.scale.factor <- sqrt(nuc.ratio.df$soma.area)
soma.scale.factor <- soma.scale.factor / max(soma.scale.factor)
soma.scale.factor2 <- nuc.ratio.df$soma.area - min(nuc.ratio.df$soma.area)
soma.scale.factor2 <- soma.scale.factor2 / max(soma.scale.factor2)
plot(nuc.ratio.df$depth_from_pia, nuc.ratio.df$nuc.vol.ratio, pch = 19,
     cex = 2 * soma.scale.factor,
     col = grey(1 - soma.scale.factor2), 
     main = "Nuc volume ratio")
lm1 <- lm(nuc.vol.ratio ~ depth_from_pia, nuc.ratio.df)
summary(lm1)
abline(lm1)

plot(nuc.ratio.df$depth_from_pia, nuc.ratio.df$nuc.area, 
     cex = 3 * nuc.ratio.df$nuc.vol.ratio, 
     main = "Nuc area")
lm1 <- lm(nuc.area ~ depth_from_pia, nuc.ratio.df)
summary(lm1)
abline(lm1)

plot(nuc.ratio.df$depth_from_pia, nuc.ratio.df$soma.area, 
     cex = 3 * nuc.ratio.df$nuc.vol.ratio, 
     main = "Soma area")
lm1 <- lm(soma.area ~ depth_from_pia, nuc.ratio.df)
summary(lm1)
abline(lm1)


g.nucprop <- ggplot(nuc.ratio.df, aes(x = depth_from_pia, y = nuc.vol.ratio)) +
  facet_wrap(~ soma.area.bin) +
  geom_point() +
  theme_bw()
plot(g.nucprop)


#### Final plot ####
g1 <- ggplot(nuc.ratio.df, aes(y = -depth_from_pia, x = nuc.ratio.df$nuc.vol.ratio,
                               size = sqrt(soma.area),
                               color = soma.area)) +
  geom_point(show.legend = FALSE) +
  xlab("Nuclear proportion") +
  ylab("Depth from pia") +
  ggtitle("Nuclear proportion of soma") +
  scale_color_gradient(low = "grey99", high = "black") +
  theme_bw()
plot(g1)
