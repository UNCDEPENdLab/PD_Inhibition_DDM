gen_model_strings <- function(diag_path){
  # browser()
  
  #this approach assumes that model traces are saved in the input diag_path and are generally labeled: "<model_string>_traces_<trace_number>.csv"
  
  all_files <- list.files(diag_path)
  all_files <- all_files[grepl("traces", all_files)] #remove DIC and GR statistics
  all_files <- sub("_traces", "", all_files) 
  all_files <- sub("_[0-9]", "", all_files) 
  all_files <- unique(sub(".csv", "", all_files))
  
}
