library("reshape")

# plot G-matrix in a constant environment
the.data <- read.table("summary_gmatrix_maternal.csv",sep=";",header=T)

# in this dataset, we have varied
# - mu_m
# - rmu
# - omega_12
# - omega_11

# restructure dataset so that I can use bwplot to plot the different G-matrices
g.dat <- melt(the.data,measure.vars=c("G11","G12","G22"))

g.dat <- g.dat[g.dat$omega_12 > 0 & g.dat$omega_11 > 0,]


for (omega_11_i in sort(unique(g.dat$omega_11)))
{
    print(paste("omega: ",omega_11_i))
    pdf(paste("g_constant",omega_11_i,".pdf",sep=""),width=25,height=12)
    print(bwplot(value ~ variable | rmu * omega_12, 
                    auto.key=T,
                    pch="|", 
                    panel=panel.superpose, 
                    groups=mu_m, 
                    data=g.dat[g.dat$omega_11 == omega_11_i,],
                    box.width=1/3, 
                    panel.groups=function(x, y, ..., group.number) {
                        panel.bwplot(x + (group.number-1.5)/3, y, ...)
                        print(length(y))
                        print(paste("next...",group.number))
                        print(y)
                    },
                    strip=function(...,which.given,factor.levels,strip.levels) { print(paste("wchi given: ", which.given));  strip.default(...,which.given,factor.levels,strip.levels=T) }
                )
        )
    dev.off()
}

# restructure dataset so that I can use bwplot to plot the different G-matrices
ev.dat <- melt(the.data,measure.vars=c("ev1","ev2"))

ev.dat <- ev.dat[ev.dat$omega_12 > 0 & ev.dat$omega_11 > 0,]


for (omega_11_i in sort(unique(ev.dat$omega_11)))
{
    print(paste("omega: ",omega_11_i))
    pdf(paste("g_constant_evals_",omega_11_i,".pdf",sep=""),width=25,height=12)
    print(bwplot(value ~ variable | rmu * omega_12, 
                    auto.key=T,
                    pch="|", 
                    panel=panel.superpose, 
                    groups=mu_m, 
                    data=ev.dat[ev.dat$omega_11 == omega_11_i,],
                    box.width=1/3, 
                    panel.groups=function(x, y, ..., group.number) {
                        panel.bwplot(x + (group.number-1.5)/3, y, ...)
                        print(length(y))
                        print(paste("next...",group.number))
                        print(y)
                    },
                    strip=function(...,which.given,factor.levels,strip.levels) { print(paste("wchi given: ", which.given));  strip.default(...,which.given,factor.levels,strip.levels=T) }
                )
        )
    dev.off()
}
