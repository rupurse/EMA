#Plot tongue tip height with red line at upper incisor sensor height
TTHeight<-ggplot(EMAnorm, aes(time_s, TTy, alpha=TTx))
TTHeight + geom_line(size=1) + geom_hline(aes(yintercept=0), colour="red") + theme_bw() +
  scale_x_continuous(limits = c(3.9, 4.37)) + scale_y_continuous(limits = c(-15, 15)) +
  ylab("") + xlab("") + theme(legend.position="none")
