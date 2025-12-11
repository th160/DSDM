library(depmixS4)

#load the trained model
final_params <- readRDS("../model/DNase-seq/Dnase_model_parameters.rds")
#finalModel <- readRDS("Dnase-seq_model.rds")

#load the extracted DNA shape features of yours DNA sequence
df3 <- readRDS("base_rs445.rds")
ntimes_trimmed <- rep(19, 2)

#Model set
mod <- depmix(list(MGW~1+MGWlag1+MGWlag2+MGWlag3,
                   EP~1+EPlag1,
                   ProT~1+ProTlag1,
                   Roll~1,
                   HelT~1+HelTlag1),
              data=df3,
              transition = ~1,
              nstates=9,
              family=list(gaussian(),
                          gaussian(),gaussian(),gaussian(),gaussian()), 
              ntimes=ntimes_trimmed)

#modNew <- setpars(mod,getpars(finalModel)) #This is for the "Dnase-seq_model.rds"
modNew <- setpars(mod,final_params)

viterbiResults = viterbi(modNew)
fb_info <- forwardbackward(modNew)

# save the prediction results
all_data <- list(
  viterbiResults = viterbiResults,
  fb_info = fb_info
)

saveRDS(all_data, "viterbiResults_base_rs445.rds")