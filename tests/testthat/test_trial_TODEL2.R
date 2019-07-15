
dat$ind_ec

for(sel in 1:107){

  fami <- dat$cluster_id[dat$ind_ec[sel, 1]:dat$ind_ec[sel, 2]]
  subj <- dat$subj_id[dat$ind_ec[sel, 1]:dat$ind_ec[sel, 2]]
  time <- dat$time[dat$ind_ec[sel, 1]:dat$ind_ec[sel, 2]]

  sigmasq0 <- 0.5
  sigmasq1 <- 2
  sigmasq2 <- 4
  delta1 <- 1
  delta2 <- 4

  out_ext1 <- out_ext2 <- out_ext3 <- matrix(0,nrepeat_ec[sel],nrepeat_ec[sel])
  for(i in 1:nrepeat_ec[sel]){
    for(j in 1:nrepeat_ec[sel]){
      if(fami[i] == fami[j]){
        out_ext1[i, j] <- sigmasq0
      }
      if(fami[i] == fami[j]){
        #if(subj[i] != subj[j]){
          out_ext2[i, j] <- sigmasq1 * exp(-abs(time[i] - time[j])/delta1)
        #}
      }
      if(subj[i] == subj[j]){
        out_ext3[i, j] <- sigmasq2 * exp(-abs(time[i] - time[j])/delta2)
      }
    }
  }

  library(matrixcalc)

  print(is.positive.definite(out_ext1))
  print(is.positive.definite(out_ext2))
  print(is.positive.definite(out_ext3))

}


