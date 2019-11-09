library("randomcoloR")
# Fig 1B
#chr.col = randomColor(30)
chr.col = distinctColorPalette(30)
chr.col1 = chr.col

# labels
dp.labs = c("1(Z)",2:30)
sl.labs = c(1:30, "31(Z)")

# Dp cordinates
ideo.dp = read.table("../# Danaus_HiRise/Dapl_Zhan_v3_HiC.chr.seq_length_dp.txt", header = F, sep = "", as.is = T)
a = 0
for (i in 1:30) {
  ideo.dp$Start[i] = a 
  ideo.dp$End[i] = ideo.dp$Start[i] + ideo.dp[i,2]
  a = ideo.dp$End[i] + 2000000
}
ideo.dp$col = chr.col

# Sl cordinates
ideo.sl = read.table("../# Danaus_HiRise/karyotype.sl.txt", header = F, sep = "", as.is = T)
sl.order = rev(c(3,13,19,10,4,26,12,24,18,8,5,11,28,1,17,23,6,25,22,7,20,14,21,27,9,15,16,29,2,31,30))
ideo.sl.r = ideo.sl[sl.order,c(3,6)]
a = 0
for (i in 1:31) {
  ideo.sl.r$Start[i] = a 
  ideo.sl.r$End[i] = ideo.sl.r$Start[i] + ideo.sl.r[i,2]
  a = ideo.sl.r$End[i] + 2000000
}

# Synteny
#ideo.sl$V6 = ideo.sl$V6*0.6
Sl.Dp = read.table("../# Danaus_HiRise/synteny.SlDp.0.8.txt", header = F, sep = "", as.is = T)
#plot(0, xlim = c(0, 325000000), ylim = c(0, 50), type = "n", yaxt = "n", xaxt = "n", xlab = "", ylab = "")

Sl.Dp_1 = subset(Sl.Dp, V4!="dp1")
Sl.Dp_Z = subset(Sl.Dp, V4=="dp1")

# re-position Dp: 
# Dp_space = (450000000-245185305-29*2000000)/2
Dp_space = 80000000

ideo.dp$Start_r = ideo.dp$Start + Dp_space
ideo.dp$End_r = ideo.dp$End + Dp_space

pdf("Fig_1B_Synteny_DpSl_r.pdf", width = 12, height = 5, pointsize = 16)
par(mar=c(0, 0, 0, 1.1))
plot(0, xlim = c(0, 450000000), ylim = c(10, 50), type = "n", axes = F, xlab = "", ylab = "")
#yaxt = "n", xaxt = "n"
segments(0+Dp_space, 40, 5685560+Dp_space, 40, lwd = 5, lend = 1, col = "olivedrab")
segments(5685561+Dp_space, 40, ideo.dp[1,7], 40, lwd = 5, lend = 1, col = "darkorange")
# ideograms:Dp
for (i in 2:30) {
  segments(ideo.dp[i,6], 40, ideo.dp[i,7], 40, lwd = 5, lend = 1, col = chr.col[i])
}
# labels
for (i in 1:30) {
  text((ideo.dp[i,6]+ideo.dp[i,7])/2, 42, dp.labs[i], col = "black", cex = 0.5)
}

# ideograms:Sl
for (i in 1:31) {
  segments(ideo.sl.r[i,3], 20, ideo.sl.r[i,4], 20, lwd = 5, lend = 1, col = "grey50")
}
for (i in 1:31) {
  text((ideo.sl.r[i,3]+ideo.sl.r[i,4])/2, 18, rownames(ideo.sl.r)[i], col = "black", cex = 0.5)
}

# synteny: autosomes
for (i in 1:nrow(Sl.Dp_1)) {
  segments(ideo.dp$Start_r[ideo.dp$V1==Sl.Dp_1[i,4]]+Sl.Dp_1[i,5], 39.5, 
           ideo.sl.r$Start[ideo.sl.r$V3==Sl.Dp_1[i,1]]+Sl.Dp_1[i,2], 20.5, lwd = 0.2, lend = 1, 
           col = ideo.dp$col[ideo.dp$V1==Sl.Dp_1[i,4]])
}
# synteny: Z
for (i in 1:nrow(Sl.Dp_Z)) {
  if (Sl.Dp_Z[i,5]<5685560) {
    b = "olivedrab"
  } else {
    b = "darkorange"
  }
  segments(ideo.dp$Start_r[ideo.dp$V1==Sl.Dp_Z[i,4]]+Sl.Dp_Z[i,5], 39.5, 
           ideo.sl.r$Start[ideo.sl.r$V3==Sl.Dp_Z[i,1]]+Sl.Dp_Z[i,2], 20.5, lwd = 0.2, lend = 1, col = b)
}

# species title
text(ideo.dp[30,4]/2 + Dp_space, 45, "Danaus plexippus", font = 3, cex = 0.8)
text(ideo.sl.r[31,4]/2, 15, "Spodoptera litura", font = 3, cex = 0.8)

dev.off()

### Scaled
# scaling factor
245185305/438940887
[1] 0.5585839
s_fac = 0.6
ideo.sl.r$Start_s[1] =  ideo.sl.r$Start[1]
ideo.sl.r$End_s[1] =  ideo.sl.r$End[1]*s_fac
for (i in 2:31) {
  ideo.sl.r$Start_s[i] =  ideo.sl.r$End_s[i-1] + 2000000
  ideo.sl.r$End_s[i] =  ideo.sl.r$Start_s[i] + ideo.sl.r[i,2]*s_fac
}

pdf("Fig_1B_Synteny_DpSl_s.pdf", width = 12, height = 5, pointsize = 16)
par(mar=c(0, 0, 0, 1.1))
plot(0, xlim = c(0, 300000000), ylim = c(10, 50), type = "n", axes = F, xlab = "", ylab = "")
segments(0, 40, 5685560, 40, lwd = 5, lend = 1, col = "olivedrab")
segments(5685561, 40, ideo.dp[1,4], 40, lwd = 5, lend = 1, col = "darkorange")
# ideograms:Dp
for (i in 2:30) {
  segments(ideo.dp[i,3], 40, ideo.dp[i,4], 40, lwd = 5, lend = 1, col = chr.col[i])
}
# labels
dp.labs = c("1(Z)",2:30)
for (i in 1:30) {
  text((ideo.dp[i,3]+ideo.dp[i,4])/2, 42, dp.labs[i], col = "black", cex = 0.5)
}

# ideograms:Sl
for (i in 1:31) {
  segments(ideo.sl.r[i,5], 20, ideo.sl.r[i,6], 20, lwd = 5, lend = 1, col = "grey50")
}
# labels
#sl.labs = c(1:30, "31(Z)")
for (i in 1:31) {
  text((ideo.sl.r[i,5]+ideo.sl.r[i,6])/2, 18, rownames(ideo.sl.r)[i], col = "black", cex = 0.5)
}

# synteny: autosomes
Sl.Dp_1 = subset(Sl.Dp, V4!="dp1")
for (i in 1:nrow(Sl.Dp_1)) {
  segments(ideo.dp$Start[ideo.dp$V1==Sl.Dp_1[i,4]]+Sl.Dp_1[i,5], 39.5, 
           ideo.sl.r$Start_s[ideo.sl.r$V3==Sl.Dp_1[i,1]]+Sl.Dp_1[i,2]*s_fac, 20.5, lwd = 0.2, lend = 1, 
           col = ideo.dp$col[ideo.dp$V1==Sl.Dp_1[i,4]])
}
# synteny: Z
Sl.Dp_Z = subset(Sl.Dp, V4=="dp1")
for (i in 1:nrow(Sl.Dp_Z)) {
  if (Sl.Dp_Z[i,5]<5685560) {
    b = "olivedrab"
  } else {
    b = "darkorange"
  }
  segments(ideo.dp$Start[ideo.dp$V1==Sl.Dp_Z[i,4]]+Sl.Dp_Z[i,5], 39.5, 
           ideo.sl.r$Start_s[ideo.sl.r$V3==Sl.Dp_Z[i,1]]+Sl.Dp_Z[i,2]*s_fac, 20.5, lwd = 0.2, lend = 1, col = b)
}

# species
text(ideo.dp[30,4]/2, 45, "Danaus plexippus", font = 3, cex = 0.8)
text(ideo.sl.r[31,6]/2, 15, "Spodoptera litura", font = 3, cex = 0.8)

dev.off()
