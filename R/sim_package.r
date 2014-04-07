###################################
## This contains 1 main function "popRangeSim"
## that is a wrapper for a python program 
## that simulates complex, stochastic demographic models
###################################

popRangeSim <- function(world, popSize, rMean = 0, rVar = 0, A=0, K=100, catProb=0, 
                          diploid=TRUE, nGens=100, migration=0, SNP_model = 1,
                          h = 0.5, s = 0, gamma_alpha = NULL, gamma_beta = NULL,
                          gSize = NULL, mutRate = NULL, nSNPs = NULL, SNPs_starting_freq = NULL,
                          GENEPOP = FALSE, GENELAND = FALSE, PLINK = FALSE, 
                          outfile = "", infile=NULL, recordTrag = 0, sDiff=NULL) {
  
  options("scipen"=100, "digits"=4)
  
  a = file.path(tempdir(), 'world.txt')
  write.table(world, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  input_call = paste(c('-p ', a, ' --outfile ', outfile, ' --recordTrag ', recordTrag,
                ' --SNP_model ', as.character(SNP_model), ' --h_input ', as.character(h),
                ' --nGens ', as.character(nGens)), collapse='')
                #input_call = paste(c("-p world.txt --outfile ", outfile, ' --recordTrag ', recordTrag,
   #                    ' --SNP_model ', as.character(SNP_model), ' --h_input ', as.character(h),
   #                    ' --nGens ', as.character(nGens)), collapse='')

  if (is.null(infile) == FALSE) { input_call = paste(input_call,' --infile ', as.character(infile)) }
  
  if (is.null(sDiff) == FALSE){
    a = file.path(tempdir(), 'sDiff.txt')
    input_call = paste(c(input_call, ' --sDiff "', a, '"'), collapse='')
    write.table(sDiff, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=',')
  }

  #SNP_model...1 == 'fixed'. 2 == 'snm'
  
  if (GENELAND == TRUE) { input_call = paste(input_call, '--GENELAND', 'TRUE') }
  else { outfile_input = paste(input_call, '--GENELAND', 'FALSE') }
  if (GENEPOP == TRUE) { input_call = paste(input_call, '--GENEPOP', 'TRUE') }
  else { outfile_input = paste(input_call, '--GENEPOP', 'FALSE') }
  if (PLINK == TRUE) { input_call = paste(input_call, '--PLINK', 'TRUE') }
  else { outfile_input = paste(input_call, '--PLINK', 'FALSE') }

  if ( is.null(gamma_alpha) & is.null(gamma_beta) ){ 
    if (is.null(nSNPs) == FALSE){
      input_call = paste(input_call, '--nSNPs_input', as.character(nSNPs), '--SNPs_starting_freq', 
                         as.character(SNPs_starting_freq))
    }
    if (class(s) == "matrix"){
      a = file.path(tempdir(), 's.txt')
      input_call = paste(c(input_call, ' --s_mat "', a, '"'), collapse='')
      write.table(s, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=',')
    }
    else if (class(s) == 'numeric') { 
      input_call = paste(input_call, '--s_input', as.character(s))
    }
  }
  else {
    input_call = paste(input_call, "--gamma_alpha_input", as.character(gamma_alpha), 
                     "--gamma_beta_input", as.character(gamma_beta))
  }
 
  if (is.null(gSize) == FALSE){
    input_call = paste(input_call, '--gSize_input', as.character(gSize), '--mutRate_input', as.character(mutRate))
  }
  
  #Migration
  if (class(migration) == "numeric") { mig_input = paste("--mig_int", as.character(migration))}
  else if (class(migration) == "matrix") {
    a = file.path(tempdir(), 'mig.txt')
    mig_input = paste(c('--mig_mat "', a, '"'), collapse='')
    write.table(migration, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  }
  
  #Ne
  if (class(popSize) == "numeric"){ Ne_input = paste('-n', popSize) }
  else if (class(popSize) == "matrix") {
    a = file.path(tempdir(), 'ne.txt')
    Ne_input = paste(c('-o "', a, '"'), collapse='') #Note this is Ne_Mat
    write.table(popSize, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  }
  else{ print("ERROR.  Ne must be either of class 'numeric' or 'matrix'") }
  
  #growth
  if (class(rMean) == "numeric"){ rMean_input=paste('--rMean_int', rMean) }
  else if (class(rMean) == "matrix") {
    a = file.path(tempdir(), 'rMean.txt')
    rMean_input = paste(c('--rMean_mat "', a, '"'), collapse='')
    write.table(rMean, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  }
  else{ print("ERROR.  rMean must be either of class 'numeric' or 'matrix' ") }

  if (class(rVar) == "numeric"){ rVar_input=paste('--rVar_int', rVar) }
  else if (class(rVar) == "matrix") {
    a = file.path(tempdir(), 'rVar.txt')
    rVar_input= paste(c('--rVar_mat "', a, '"'), collapse='')
    write.table(rVar, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  }
  else{ print("ERROR.  rVar must be either of class 'numeric' or 'matrix' ") }
  
  if (class(A) == "numeric"){ A_input = paste('--A_int', A) }
  else if (class(A) == "matrix") {
      a = file.path(tempdir(), 'A.txt')
      A_input = paste(c('--A_mat "', a, '"'), collapse='')
      write.table(A, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  }
  else{ print("ERROR. A must be either of class 'numeric' or 'matrix' ")}
  
  if (class(K) == "numeric"){ K_input = paste('--K_int', K) }
  else if (class(K) == "matrix") {
    a = file.path(tempdir(), 'K.txt')
    K_input = paste(c('--K_mat "', a, '"'), collapse='')
    write.table(K, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  }
  else{ print("ERROR. K must be either of class 'numeric' or 'matrix' ") }
  
  if (class(catProb) == "numeric"){ catProb_input = paste('--catProb_int', catProb)}
  else if (class(catProb) == "matrix") {
    a = file.path(tempdir(), 'catProb.txt')
    catProb_input = paste(c('--catProb_mat "', a, '"'), collapse='')
    write.table(catProb, a, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
  }
  else{ print("ERROR. catProb must be either of class 'numeric' or 'matrix' ")}
  
  if (diploid == TRUE){ diploid_input = paste('--diploid', 'True')}
  else { diploid_input = paste('--diploid', 'False')}

  
  #In my case, input 1 will be the matrix of the populations
  #command <- paste("python FS_3.py ", world_input, Ne_input, rMean_input, rVar_input,
  #                 A_input, K_input, catProb_input, diploid_input, nGens_input, mig_input,
  #                 SNP_model_input, h_input, s_input, outfile_input, bb_input, recordTrag_input)
  fileLoc = paste(c('"', system.file("popRange_main.py", package="popRange"), '"'), collapse='')
  command <- paste("python", fileLoc, input_call, Ne_input, rMean_input, rVar_input,
                   A_input, K_input, catProb_input, diploid_input, mig_input)
  #command <- paste("python popRange/inst/popRange_main.py", input_call, Ne_input, rMean_input, rVar_input,
  #                 A_input, K_input, catProb_input, diploid_input, mig_input)
  system(command)
}
