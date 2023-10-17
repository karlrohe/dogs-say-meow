# install.packages("tidytext")
library(tidytext)
bff_with_unigrams = function(loadings, vector_of_strings, num_best=10, cut_off=10){
  #  run vsp::bff (best feature function) on the 
  #  document term matrix created from vector_of_strings.
  
  #  ESSENTIAL: the rows of loadings must perfectly align with the elements of vector_of_strings.
  
  
  # make vector_of_strings a document-term matrix (dtm_sparse)
  vector_of_strings = paste("xyzJunkGetsDiscardedButFixesNullStrings", vector_of_strings)
  text_data <- tibble(doc_id = 1:length(vector_of_strings), text = vector_of_strings)
  
  tidy_text <- text_data %>%
    unnest_tokens(word, text) %>%
    anti_join(stop_words)
  
  dtm_sparse <- tidy_text %>%
    count(doc_id, word) %>%
    cast_sparse(doc_id, word, n)
  
  drop_these= which(colSums(dtm_sparse)<cut_off)
  dtm_sparse = dtm_sparse[,-drop_these]
  
  num_cols <- ncol(loadings)
  max_digits <- nchar(num_cols) # Determines the number of digits in the largest column number
  
  colnames(loadings) <- sprintf("factor_%0*d", max_digits, 1:num_cols)
  
  
  out= vsp::bff(loadings, dtm_sparse, num_best)
  return(out)
  
}