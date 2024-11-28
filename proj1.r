#╔════════════════════════════════════════════════╗
#║     MATH111762 - Group 14 - Assignment 1       ║
#╠════════════════════╦════════════╦══════════════╣
#║        Name        ║ Student ID ║ Contribution ║
#╠════════════════════╬════════════╬══════════════╣
#║ Alexander Cheetham ║ s2638244   ║ 1/3          ║
#║ Thalia Andreou     ║ s2665822   ║ 1/3          ║
#║ Anni Tziakouri     ║ s2704516   ║ 1/3          ║
#╚════════════════════╩════════════╩══════════════╝ 

# Team member Contributions:
# q3: Anni Tziakouri 
# q4: Thalia Andreou 
# q5: Alexander Cheetham 
# q6: Anni Tziakouri 
# q7: Thalia Andreou 
# q8: Alexander Cheetham
# q9: Thalia Andreou 
# q10: Alexander Cheetham
# General Commenting: Everyone
# Test Cases: Alexander Cheetham
# Final checks : Anni Tziakouri


# This R script implements a pth-order Markov model to generate random text 
# based on word patterns from James Joyce’s Ulysses. 
# The pth-order Markov mode works by counting the 
# number of occurrences for a given set of "token"s in a row in that text,
# and using these frequencies as probabilities to predict the token after.

# For example with the text "my name is John, my name is James, my name is Jacob"
# a markov model of order 3 conditioned on the tokens "my name is" generates
# names John, James and Jacob with probability 1/3 each.



# set seed for reproducibility
# set.seed(12345); 


# Question 3: Load in the Ulysses text file.
# read the file into R between specified lines
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8")
# replaces all occurrences of the substring "_(" with an empty string
a <- gsub("_(","",a,fixed=TRUE) 

# Question 4: Preprocessing, separate punctuation into its own element on new vector
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

# Question 5: remove punctuation

# use punctuation as defined in rubric
punctuations <- c(",",".", ";", "!", ":", "?")

# copy a to perform removing of punctuation (avoid changing the original vector) 
a_removed_punct <-a

# for each punctuation mark remove it and assign as a new element after word 
# it was attached to
for (punct in punctuations){
  a_removed_punct<-split_punct(vec=a_removed_punct,punct=punct)}
  

# Q6: create a vector of the m most common unique words

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


# Q6 - extra graph to visualize difference with adjusting the 
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


#Question 7 : create matrix M
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

#set mlag as 4
mlag<-4

# instantiate index vector, where indexes of b are matched with lowercase_a
index_vector <- match(lowercase_a,b)

#apply function to create M
M<-create_common_token_sequences(mlag,index_vector)
#print(M)


#Question 8: simulate nw-word sections from your model

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
      if (length(unique(row_matches[,j+1]))==1|
          (length(unique(row_matches[,j+1]))==2&
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
      sampling_vector <- unique(sampling_vector[!is.na(sampling_vector)])
      # check if prob argument provided
      if(!is.null(prob)){
        # create a subset of weights for sample based on unique tokens
        # in sampling vector.
        subset_prob<-prob[unique(sampling_vector)]
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
# q8 actually generate some simulated text.

# inference the markov model 
simulated_text <- inference_model(50,index_vector,M=M,prob=NULL)
# print the simulated text with cat
print("------------START N-gram Simulated sentence----------------")
cat(b[simulated_text], sep=" ")
print("------------END N-gram Simulated sentence----------------")

# q9: compare to 50 randomly sampled words
corresponding_freqs <- frequencies[largest_indices]
random_samples <- sample(b, size = 50, replace = TRUE,prob=corresponding_freqs)
print("------------START Randomly Sampled Simulated sentence----------------")
cat(random_samples, " ")
print("------------END Randomly Sampled Simulated sentence----------------")


#q10: format the simulated text so commonly capitalised words are capitalised
#       and punctuation is printed in expected positions.

# capitalise all first letters in b
capitalise_first_letter <- function(string) {
  paste0(toupper(substring(string, 1, 1)), tolower(substring(string, 2)))
}
b_capitalized<-sapply(b, capitalise_first_letter)

# find the matches for capitalised versions of the words in B
b_cap_match <- match(a_removed_punct,b_capitalized)
b_cap_match <- b_cap_match[!is.na(b_cap_match)]
b_cap_freq <- tabulate(b_cap_match,nbins=length(b_capitalized))

# calculate the ratio of capitalised to non-capitalised versions of a given word
# using the corresponding freqs (totals)
freq_cap_ratio <- b_cap_freq / corresponding_freqs 

# get the indexes of the words that are more likely (prob 0.5+)
# or even chances to be capitalised
indexes_capitals <- which(freq_cap_ratio >= 0.5)
b_modified <- b
# capitalise all words likely to be capitalised in common word list
b_modified[indexes_capitals] <- tools::toTitleCase(b)[indexes_capitals]

write_sentence <- function(words) {
  # Function: write_sentence
  # Description: This function takes a vector of words and punctuation and 
  #             prints a sentence with correct punctuation placement.
  #              
  # Inputs:
  #   - words: vector of words and punctuation
  # Outputs:
  #   - print to console sentence with correct punctuation
  sentence <- ""
  for (i in 1:length(words)) {
    if (words[i] %in% punctuations) {
      # If the element is a punctuation mark, add it directly without space
      sentence <- paste0(sentence, words[i], " ")
    } else {
      # Add a space before the word if it's not the first word
      if (i > 1 && !words[i-1] %in% punctuations) {
        sentence <- paste0(sentence, " ")
      }
      # Add the word
      sentence <- paste0(sentence, words[i])
    }
  }
  
  # Use cat() to print the sentence without quotes
  cat(sentence)
}

print("------------START Capitalised N-gram Simulated sentence----------------")
write_sentence(b_modified[simulated_text])
print("------------END Capitalised N-gram Simulated sentence----------------")

print("------------START Capitalised Random Sampled Simulated sentence----------------")
cat(write_sentence(random_samples), " ")
print("------------END Capitalised Randomly Sampled Simulated sentence----------------")


