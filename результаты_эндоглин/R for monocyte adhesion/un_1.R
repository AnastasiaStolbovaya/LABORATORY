m_means$mab <- factor(m_means$mab, levels=c( "control", "TGF-beta", "3A7+TGF-beta" ,  "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta"))

levels(m_means$mab) <- c( "К", "TGF-beta", "Ig1+TGF-beta" , "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta")



fig1_count_U937_TGF_sum <- ggplot(m_means, aes(x = mab, y = emmean_resp)) + 
        geom_col(fill = "gray", color = "black") + 
        geom_errorbar(aes(ymin = asymp.LCL_resp, ymax = asymp.UCL_resp), width = 0.2)+
        theme_bw()  +
        theme(axis.text.x = element_text(angle = 30, hjust=1))+
        ylab("")+
        xlab("")+
        geom_segment(aes(x=c("Ig1+TGF-beta"), xend=c("4C9+TGF-beta"), y=1200, yend=1200))+
        geom_segment(aes(x=c("Ig1+TGF-beta"), xend=c("4E4+TGF-beta"), y=1210, yend=1210))+
        geom_segment(aes(x=c("Ig1+TGF-beta"), xend=c("5H7+TGF-beta"), y=1220, yend=1220))+
        geom_text(aes(x = "2C8+TGF-beta", y = 1230, label = "***"), size = 7)

