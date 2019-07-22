# After running 'waveqtl_hmt_test_calc_sumlog_init_val.R
# 20 iterations max for each case should be sufficient

### seq(0.2, 0.8, 0.2) ----

plot(logl_list[[1]][-1], col = "black",type = "l", ylab = "LogL", xlab = "# iterations")
lines(logl_list[[2]][-1], col = "blue")
lines(logl_list[[3]][-1], col = "red")
lines(logl_list[[4]][-1], col = "green")
lines(logl_list[[5]][-1], col = "black", lty = 2)
lines(logl_list[[6]][-1], col = "blue", lty = 2)
lines(logl_list[[7]][-1], col = "red", lty = 2)
lines(logl_list[[8]][-1], col = "green", lty = 2)
title("0.7/0.3 bad case 1")
legend("bottomright",pt.cex = 0.8
       ,legend = c("0.2/0.2"
                   ,"0.4/0.4"
                   ,"0.6/0.6"
                   ,"0.8/0.8"
                   ,"0.2/0.8"
                   ,"0.4/0.6"
                   ,"0.6/0.4"
                   ,"0.8/0.2")
       ,col = rep(c("black","blue","red","green"),2)
       ,lty = c(rep(1,4),rep(2,4)))

# Level 2
lvl2 <- lapply(eps_vals, function(x){
  z <- rbindlist(lapply(x, function(y){
    y[2,]
  }))
})
# Level 3
lvl3 <- lapply(eps_vals, function(x){
  z <- rbindlist(lapply(x, function(y){
    y[3,]
  }))
})
for(i in 1:length(lvl2)){
  lvl2[[i]][,"paramSet" := i]
  lvl2[[i]][,"iteration" := 1:.N]
  lvl2[[i]][,"logL" := logl_list[[i]][-1]]
  lvl3[[i]][,"paramSet" := i]
  lvl3[[i]][,"iteration" := 1:.N]
  lvl3[[i]][,"logL" := logl_list[[i]][-1]]
}
lvl2 <- copy(rbindlist(lvl2))
lvl3 <- copy(rbindlist(lvl3))

lvl2[,"paramSet" := factor(paramSet
                           ,levels = 1:8
                           ,labels = c("0.2/0.2"
                                       ,"0.4/0.4"
                                       ,"0.6/0.6"
                                       ,"0.8/0.8"
                                       ,"0.2/0.8"
                                       ,"0.4/0.6"
                                       ,"0.6/0.4"
                                       ,"0.8/0.2"))]
lvl3[,"paramSet" := factor(paramSet
                           ,levels = 1:8
                           ,labels = c("0.2/0.2"
                                       ,"0.4/0.4"
                                       ,"0.6/0.6"
                                       ,"0.8/0.8"
                                       ,"0.2/0.8"
                                       ,"0.4/0.6"
                                       ,"0.6/0.4"
                                       ,"0.8/0.2"))]

library(ggplot2)
ggplot(lvl2[,.(`11`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `11`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - tying lvl 2") +
  scale_color_manual(values = rep(c("black","blue","red","green"),2)) +
  scale_linetype_manual(values = c(rep("solid",4),rep("dashed",4))) +
  xlim(c(0,20))
ggplot(lvl2[,.(`10`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `10`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - tying lvl 2") +
  scale_color_manual(values = rep(c("black","blue","red","green"),2)) +
  scale_linetype_manual(values = c(rep("solid",4),rep("dashed",4))) +
  xlim(c(0,20))
ggplot(lvl3[,.(`11`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `11`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - tying lvl 3") +
  scale_color_manual(values = rep(c("black","blue","red","green"),2)) +
  scale_linetype_manual(values = c(rep("solid",4),rep("dashed",4)))
ggplot(lvl3[,.(`10`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `10`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - tying lvl 3") +
  scale_color_manual(values = rep(c("black","blue","red","green"),2)) +
  scale_linetype_manual(values = c(rep("solid",4),rep("dashed",4)))

### seq(0.6, 0.8, 0.1) ----
# Bottom level initialised for pr 0 transition to 1 state (as per data)

plot.new()
png("analysis/sim2_debug/s07_03_logL_initVals_bottom0.png", width = 900, height = 600, pointsize = 20)
plot(logl_list[[1]][-1], col = "black",type = "l", ylab = "LogL", xlab = "# iterations")
lines(logl_list[[2]][-1], col = "blue")
lines(logl_list[[3]][-1], col = "red")
lines(logl_list[[4]][-1], col = "black", lty = 2)
lines(logl_list[[5]][-1], col = "blue", lty = 2)
lines(logl_list[[6]][-1], col = "red", lty = 2)
title("0.7/0.3 bad case 1 - bottom level init 0")
legend("bottomright",pt.cex = 0.8
       ,legend = c("0.6/0.6"
                   ,"0.7/0.7"
                   ,"0.8/0.8"
                   ,"0.6/0.4"
                   ,"0.7/0.3"
                   ,"0.8/0.2")
       ,col = rep(c("black","blue","red"),2)
       ,lty = c(rep(1,3),rep(2,3)))
dev.off()


# Level 2
lvl2 <- lapply(eps_vals, function(x){
  z <- rbindlist(lapply(x, function(y){
    y[2,]
  }))
})
# Level 3
lvl3 <- lapply(eps_vals, function(x){
  z <- rbindlist(lapply(x, function(y){
    y[3,]
  }))
})
for(i in 1:length(lvl2)){
  lvl2[[i]][,"paramSet" := i]
  lvl2[[i]][,"iteration" := 1:.N]
  lvl2[[i]][,"logL" := logl_list[[i]][-1]]
  lvl3[[i]][,"paramSet" := i]
  lvl3[[i]][,"iteration" := 1:.N]
  lvl3[[i]][,"logL" := logl_list[[i]][-1]]
}
lvl2 <- copy(rbindlist(lvl2))
lvl3 <- copy(rbindlist(lvl3))

lvl2[,"paramSet" := factor(paramSet
                           ,levels = 1:6
                           ,labels = c("0.6/0.6"
                                       ,"0.7/0.7"
                                       ,"0.8/0.8"
                                       ,"0.6/0.4"
                                       ,"0.7/0.3"
                                       ,"0.8/0.2"))]
lvl3[,"paramSet" := factor(paramSet
                           ,levels = 1:6
                           ,labels = c("0.6/0.6"
                                       ,"0.7/0.7"
                                       ,"0.8/0.8"
                                       ,"0.6/0.4"
                                       ,"0.7/0.3"
                                       ,"0.8/0.2"))]

library(ggplot2)

ggplot(lvl2[,.(`11`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `11`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_11 - bottom level init 0 - tying lvl 2") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3))) +
  xlim(c(0,20))
ggsave("analysis/sim2_debug/s07_03_eps11_lvl2_initVals_bottom0.png", width = 15, height = 7.5, units = "cm")

ggplot(lvl2[,.(`10`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `10`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_10 - bottom level init 0 - tying lvl 2") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3))) +
  xlim(c(0,20))
ggsave("analysis/sim2_debug/s07_03_eps10_lvl2_initVals_bottom0.png", width = 15, height = 7.5, units = "cm")

ggplot(lvl3[,.(`11`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `11`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_11 - bottom level init 0 - tying lvl 3") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3)))
ggsave("analysis/sim2_debug/s07_03_eps11_lvl3_initVals_bottom0.png", width = 15, height = 7.5, units = "cm")

ggplot(lvl3[,.(`10`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `10`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_10 - bottom level init 0 - tying lvl 3") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3)))
ggsave("analysis/sim2_debug/s07_03_eps10_lvl3_initVals_bottom0.png", width = 15, height = 7.5, units = "cm")


# Bottom level initialised at same as other init vals (not accounting for the 0 transition)

plot.new()
png("analysis/sim2_debug/s07_03_logL_initVals.png", width = 900, height = 600, pointsize = 20)
plot(logl_list[[1]][-1], col = "black",type = "l", ylab = "LogL", xlab = "# iterations")
lines(logl_list[[2]][-1], col = "blue")
lines(logl_list[[3]][-1], col = "red")
lines(logl_list[[4]][-1], col = "black", lty = 2)
lines(logl_list[[5]][-1], col = "blue", lty = 2)
lines(logl_list[[6]][-1], col = "red", lty = 2)
title("0.7/0.3 bad case 1")
legend("bottomright",pt.cex = 0.8
       ,legend = c("0.6/0.6"
                   ,"0.7/0.7"
                   ,"0.8/0.8"
                   ,"0.6/0.4"
                   ,"0.7/0.3"
                   ,"0.8/0.2")
       ,col = rep(c("black","blue","red"),2)
       ,lty = c(rep(1,3),rep(2,3)))
dev.off()


# Level 2
lvl2 <- lapply(eps_vals, function(x){
  z <- rbindlist(lapply(x, function(y){
    y[2,]
  }))
})
# Level 3
lvl3 <- lapply(eps_vals, function(x){
  z <- rbindlist(lapply(x, function(y){
    y[3,]
  }))
})
for(i in 1:length(lvl2)){
  lvl2[[i]][,"paramSet" := i]
  lvl2[[i]][,"iteration" := 1:.N]
  lvl2[[i]][,"logL" := logl_list[[i]][-1]]
  lvl3[[i]][,"paramSet" := i]
  lvl3[[i]][,"iteration" := 1:.N]
  lvl3[[i]][,"logL" := logl_list[[i]][-1]]
}
lvl2 <- copy(rbindlist(lvl2))
lvl3 <- copy(rbindlist(lvl3))

lvl2[,"paramSet" := factor(paramSet
                           ,levels = 1:6
                           ,labels = c("0.6/0.6"
                                       ,"0.7/0.7"
                                       ,"0.8/0.8"
                                       ,"0.6/0.4"
                                       ,"0.7/0.3"
                                       ,"0.8/0.2"))]
lvl3[,"paramSet" := factor(paramSet
                           ,levels = 1:6
                           ,labels = c("0.6/0.6"
                                       ,"0.7/0.7"
                                       ,"0.8/0.8"
                                       ,"0.6/0.4"
                                       ,"0.7/0.3"
                                       ,"0.8/0.2"))]

library(ggplot2)

ggplot(lvl2[,.(`11`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `11`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_11 - tying lvl 2") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3))) +
  xlim(c(0,20))
ggsave("analysis/sim2_debug/s07_03_eps11_lvl2_initVals.png", width = 15, height = 7.5, units = "cm")

ggplot(lvl2[,.(`10`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `10`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_10 - tying lvl 2") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3))) +
  xlim(c(0,20))
ggsave("analysis/sim2_debug/s07_03_eps10_lvl2_initVals.png", width = 15, height = 7.5, units = "cm")

ggplot(lvl3[,.(`11`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `11`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_11 - tying lvl 3") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3)))
ggsave("analysis/sim2_debug/s07_03_eps11_lvl3_initVals.png", width = 15, height = 7.5, units = "cm")

ggplot(lvl3[,.(`10`,paramSet,iteration,logL)]) +
  geom_line(aes(x = iteration, y = `10`, group = paramSet, colour = factor(paramSet), linetype = factor(paramSet))) +
  ggtitle("0.7/0.3 bad case 1 - eps_10 - tying lvl 3") +
  scale_color_manual(values = rep(c("black","blue","red"),2)) +
  scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3)))
ggsave("analysis/sim2_debug/s07_03_eps10_lvl3_initVals.png", width = 15, height = 7.5, units = "cm")
