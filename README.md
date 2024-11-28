## Project 1 : pth-order Markov for Text Generation

#### Overview
File name `markov_model.r` <br>
This repository contains an implementation of a text generation model based on a pth-order Markov model, using a source text. The goal of this assignment is to investigate how easily random text can be distinguished from real text by generating sequences of words based on the word patterns in the source text. The model generates random text by using word sequences from the original text, and it simulates word generation by following the statistical properties of word sequences.

#### Implementation
Pre-processing the text
Data Preparation : token frequency analysis and create mapping for efficient storage and processing
Building the Markov Model : matrix of word sequences and their associated probabilities and handle edge cases (falling back on lower-order models)
Simulate the generation of sections of text using the Markov model and compare the generated text to another model based on frequency

## Project 2 :  Covid spread modeling

#### Overview
File name `covid_analysis.r` <br>
The goal is to infer the daily incidence of Covid infections from hospital death data. The approach is based on a simulation that iteratively adjusts the guessed infection dates of Covid victims to match the observed death distribution.

#### Implementation
Model the infection-to-Death Distribution:
Initialize Infection Dates
Simulate a loop to update infection dates, compute pearson statistic and adjust infection dates to minimize P
Plot estimated incidence over time, compare it to real deaths, and mark the UK lockdown start date.


