source("simulation.R")
filename <- "output.out"
#cat("This is the output file", file=filename, append=FALSE, sep="\n")
for(i in 1:5) {
b <- simulation()
for(i in 1:4) {
   cat(b$out1[i], file=filename, append=TRUE, sep="\n")
}
for(i in 1:4) {
   cat(b$out2[i], file=filename, append=TRUE, sep="\n")
}
for(i in 1:4) {
   cat(b$pcs[i], file=filename, append=TRUE, sep=",")
}
cat("\n\n", file=filename, append=TRUE)

#for(i in 1:10) {
#   cat(simulation(), file=filename, append=TRUE, sep="\n")
#}
}
