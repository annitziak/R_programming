# This R script implements a pth-order Markov model to generate random text 
# The pth-order Markov mode works by counting the 
# number of occurrences for a given set of "token"s in a row in that text,
# and using these frequencies as probabilities to predict the token after.

# For example with the text "my name is John, my name is James, my name is Jacob"
# a markov model of order 3 conditioned on the tokens "my name is" generates
# names John, James and Jacob with probability 1/3 each.


# set seed for reproducibility
# set.seed(12345); 


# Load in text file.
a <- scan("input.txt",what="character",
          fileEncoding="UTF-8")
# replaces all occurrences of the substring "_(" with an empty string
a <- gsub("_(","",a,fixed=TRUE) 

# Preprocessing, separate punctuation into its own element on new vector
# e.g "a,b.c?" --> "a", ",", "b", ".", "c", "?"

split_punct<-function(vec,punct){
  
  # Function: split_punct
  # Description: This function splits a vector of strings based on a 
  #              given character. It inserts that character 
  #              after the word it came from into a separate element in the vector, 
  #              while removing it from the original words.
  #              
  # Inputs:
  #   - vec: A vector of strings.
  #   - punct: A character to split the vector on.
  # Outputs:
  #   - A vector where occurrences of punct have been separated
  #     into their own elements.
  
  #check variable types
  if (!is.vector(vec)) {
    stop("Error: 'vec' must be of type 'vector'")
  }
  if (!is.character(punct)) {
    stop("Error: 'punct' must be of type 'character'")
  }
  # get indexes of words with punctuation in
  # Use fixed = TRUE for literal string matching (avoid regex)
  index_punct <-grep(punct,vec,fixed=TRUE)
  print
  #if no matches, skip
  if (length(index_punct)==0){
    return(vec)
  }
  
  #remove occurrences of punctuation in input vector
  words_vec <- gsub(punct, "",vec, fixed=TRUE)
  
  # new character vector of length original + punctuation occurrences
  replacement_vector <-rep("", length(words_vec)+length(index_punct))
  
  # calculate indices for new punctuation elements to be inserted after each punctuated word
  iis <- index_punct + 1:length(index_punct)
  
  #insert words and punctuation into new vector
  replacement_vector[iis]<-paste(punct)
  replacement_vector[-iis]<-paste(words_vec)
  
  return(replacement_vector)
  
}

# remove punctuation

# use punctuation as defined in rubric
punctuations <- c(",",".", ";", "!", ":", "?")

# copy a to perform removing of punctuation (avoid changing the original vector) 
a_removed_punct <-a

# for each punctuation mark remove it and assign as a new element after word 
# it was attached to
for (punct in punctuations){
  a_removed_punct<-split_punct(vec=a_removed_punct,punct=punct)}
  

# create a vector of the m most common unique words

# convert to lowercase, 
# collect all the unique words in from vector
# and find their matched indexes to lowercase_a
lowercase_a <-tolower(a_removed_punct)
unique_words <-unique(lowercase_a)
matched_words <- match(lowercase_a,unique_words)

#calculate the frequency of each unique word
frequencies <-tabulate(matched_words)

# set m most common words and find the threshold value 
m<-1001
threshold_value<- sort(frequencies, decreasing = TRUE)[m]

# select only the indexes of the words with frequency above the threshold value
largest_indices <- which(frequencies >= threshold_value)

# instantiate b, the list of the m most common words.
b<-unique_words[largest_indices]


# extra graph to visualize difference with adjusting the 
# m most common words

# Initialize values of m and corresponding threshold values
m_values <- seq(10, 1200, by = 10) # Adjust the range as necessary
threshold_values <- numeric(length(m_values))

# Calculate threshold values for each m
for (i in seq_along(m_values)) {
  m <- m_values[i]
  threshold_values[i] <- sort(frequencies, decreasing = TRUE)[m]
}

# Plotting the thresholds and m
plot(m_values, threshold_values, type = "l", col = "blue", 
     xlab = "Number of Most Common Words (m)", 
     ylab = "Frequency - Threshold Value", 
     main = "Threshold Values vs. Number of Most Common Words (m)",
     lwd = 2)

# Findings : it can be seen that for a small number of common words the threshold 
# value of the occurrences is very high but it rapidly decreases the more words
# we decide to consider, before it plateaus after about 300 words, indicating that
# including more words does not significantly affect the minimum frequency required.


# create matrix M
# M is of size (n − p )x(p+1) and contains chronological contiguous p+1 word 
# snippets of the tokenized version of Ulysses as its rows.

create_common_token_sequences<-function(mlag,word_vec){
  
  # Function: create_common_token_sequences
  # Description: This function creates a matrix M, where each column is the 
  # word_vector and the next column contains the vector shifted down once.
  #              
  # Inputs:
  #   - mlag: integer lag value
  #   - word_vec: vector of common tokens
  # Outputs:
  #   - Matrix M: lagged matrix with shift 1
  
  n<-length(word_vec)
  #initiate new matrix with dimensions (n-mlag)x(mlag+1)
  M<-matrix(nrow=(n-mlag),ncol=(mlag+1))
  # Loop through each column index from 1 to (mlag + 1)
  for (i in 1:(mlag+1)) {
    # Extract the current column vector from word_vec with appropriate indexing
    # and assign it to matrix M
    column_i = word_vec[i:(length(word_vec) - (mlag-(i-1)))]
    M[,i] <-column_i
  }
  return(M)
}

#set mlag as 4 for now
mlag<-4

# instantiate index vector, where indexes of b are matched with lowercase_a
index_vector <- match(lowercase_a,b)

#apply function to create M
M<-create_common_token_sequences(mlag,index_vector)
#print(M)


#simulate nw-word sections from your model

inference_model<-function(nw,index_vector,M,prob=NULL){
  # Function: inference_model
  # Description: this function implements the inference of an order p markov 
  # returning a vector of simulated words             
  # Inputs:
  #   - Matrix M: As defined in q7: (n − p )x(p+1)  matrix that contains
  #               chronological contiguous p+1 word snippets as rows
  #   - word_vec: vector of m unique most common tokens
  #   - prob: a vector of probabilities that should be assigned to each word
  #           in the first column of M/index vector
  # Outputs:
  #     - final_sequence: vector of simulated words.
  
  
  
  # remove NAs from index vector and sample the first token randomly
  index_vector_non_NA <- na.omit(index_vector)
  first_token <- sample(index_vector_non_NA,size=1)
  
  #setup output vector and fill with first token
  final_sequence <-rep(NA,nw)
  final_sequence[1] <-first_token
  
  
  for (i in 2:nw){
    for (j in mlag:1) if (i>j){
      
      # print indexes 
      print(paste("i=",i,"j=",j))
      
      # generate the current ngram from which a new token should be generated
      cur_ngram <- final_sequence[(i-j):(i-1)]
      
      # create the matrix containing just the first x columns of M
      # where x is the length of cur_ngram
      get_columns = matrix(M[,1:length(cur_ngram)],
                           nrow=(length(index_vector)-mlag),
                           ncol=length(cur_ngram),
                           byrow=FALSE)
      
      # select only those rows from get_columns where the first x columns match
      # the current ngram
      
      # 1. find all elements where ngram[i]==M[i], matrix vector comparison
      # is performed columnwise in R hence the transposition
      # transpose matrix to make matrix vector comparison work.
      # converting NAs to FALSE for logical comparisons to not return NAs
      equality <- t(t(get_columns)==cur_ngram)
      equality[is.na(equality)] <- FALSE
      
      # 2. if the sum of the row == length of ngram we have every 
      # character of that ngram matched in the row and therefore we should
      # keep that row
      row_sums <- rowSums(equality)
      matched_indexes <- row_sums== length(cur_ngram)
      
      # retrieve those matched rows
      row_matches <- M[matched_indexes,]
  
      # check to see if matched rows have zero rows of matches, one row of matches
      # or all NAs. If yes, use the p-1th order model and sample randomly
      # from first column if p-1 == 0 (i.e 0th order model)
      if (!is.matrix(row_matches)) {
        # one ngram match that is a valid token
        #print("one option")
        
        # p=0 so sample randomly from index vector
        if(j==1){
          final_sequence[i]<-sample(index_vector_non_NA,size=1)
        }
        next
      } else if (nrow(row_matches) == 0) {
        # no ngram matches at all
        #print("no options")
        
        # p=0 so sample randomly from index vector
        if(j==1){
          final_sequence[i]<-sample(index_vector_non_NA,size=1)
        }
        next
      } 
      
      # check that options are not of one type
      # including edge case of token,NA which also counts as one type
      if (length((row_matches[,j+1]))==1|
          (length((row_matches[,j+1]))==2&
           sum(is.na(row_matches[,j+1]))>=1)) {
        
        # all ngram match options are one token or {token,NA} which is 
        #the same case as one token
        
        #print("All options are one type")
        
        # p=0 so sample randomly from index vector
        if(j==1){
          final_sequence[i]<-sample(index_vector_non_NA,size=1)
        }
        next
      }
      
      # select the next valid  token sampling from the j+1-th column of matches
      sampling_vector <- row_matches[,j+1]
      sampling_vector <- sampling_vector[!is.na(sampling_vector)]
      # check if prob argument provided
      if(!is.null(prob)){
        # create a subset of weights for sample based on unique tokens
        # in sampling vector.
        subset_prob<-prob[sampling_vector]
        next_token <- sample(sampling_vector, size=1,
                             prob=subset_prob,replace = TRUE)
      } 
      else {
        next_token <- sample(sampling_vector, size=1,prob=prob,replace=TRUE)
      }
      
      final_sequence[i] <- next_token 
      # print(paste("WE MADE A TOKEN",next_token))
      
      # for debugging, note when a match of greater than three characters 
      # in a row occurs.
      if (j>=3){
        print("BIGGGGG MATCH")
      }
      break
      }
  }
  return(final_sequence)
}
# now ctually generate some simulated text.

# inference the markov model 
simulated_text <- inference_model(50,index_vector,M=M,prob=NULL)
# print the simulated text with cat
print("------------START N-gram Simulated sentence----------------")
cat(b[simulated_text], sep=" ")
print("------------END N-gram Simulated sentence----------------")

# compare to 50 randomly sampled words
corresponding_freqs <- frequencies[largest_indices]
random_samples <- sample(b, size = 50, replace = TRUE,prob=corresponding_freqs)
print("------------START Randomly Sampled Simulated sentence----------------")
cat(random_samples, " ")
print("------------END Randomly Sampled Simulated sentence----------------")



