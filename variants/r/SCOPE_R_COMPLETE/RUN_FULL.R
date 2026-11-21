source("SCOPE_Exp123_FINAL.R")
# For publication: choose separate experiment calls to simplify checkpointing.
SCOPE_Exp123_FINAL("full","Exp1",c(512L,1024L,2048L),c(3,5,7),1000L)
SCOPE_Exp123_FINAL("full","Exp2",c(512L,1024L,2048L),c(3,5,7),1000L)
SCOPE_Exp123_FINAL("full","Exp3",c(512L,1024L,2048L),c(3,5,7),1000L)
