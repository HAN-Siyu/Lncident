findPotentialStartsAndStops <- function(sequence) {
        codons<- c("atg", "taa", "tag", "tga")
        for (i in 1:4) {
                codon <- codons[i]
                codonpositions <- unlist(gregexpr(codon, sequence, ignore.case = TRUE))
                if(codonpositions[1] == -1) codonpositions <- integer()
                numoccurrences <- length(codonpositions)
                if (i == 1) {
                        positions <- codonpositions
                        types <- rep(codon, numoccurrences)
                } else {
                        positions <- append(positions, codonpositions,
                                            after = length(positions))
                        types <- append(types, rep(codon, numoccurrences),
                                        after = length(types))
                }
        }
        indices <- order(positions)
        positions <- positions[indices]
        types <- types[indices]
        mylist <- list(positions, types)
        return(mylist)
}
pos_orf <- function(sequence) {
        mylist <- findPotentialStartsAndStops(sequence)
        positions <- mylist[[1]]
        types <- mylist[[2]]
        orfstarts <- numeric()
        orfstops <- numeric()
        orflengths <- numeric()
        numpositions <- length(positions)
        if (numpositions >= 2) {
                for (i in 1:(numpositions-1)) {
                        posi <- positions[i]
                        typei <- types[i]
                        found <- 0
                        while (found == 0) {
                                for (j in (i+1):numpositions) {
                                        posj  <- positions[j]
                                        typej <- types[j]
                                        posdiff <- posj - posi
                                        posdiffmod3 <- posdiff %% 3
                                        orflength <- posj - posi + 3
                                        if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0) {
                                                numorfs <- length(orfstops)
                                                usedstop <- -1
                                                if (numorfs > 0) {
                                                        for (k in 1:numorfs) {
                                                                orfstopk <- orfstops[k]
                                                                if (orfstopk == (posj + 2)) { usedstop <- 1 }
                                                        }
                                                }
                                                if (usedstop == -1) {
                                                        orfstarts <- append(orfstarts, posi, after = length(orfstarts))
                                                        orfstops <- append(orfstops, posj + 2, after = length(orfstops))
                                                        orflengths <- append(orflengths, orflength, after = length(orflengths))
                                                }
                                                found <- 1
                                                break
                                        }
                                        if (j == numpositions) { found <- 1 }
                                }
                        }
                }
        }
        indices <- order(orfstarts)
        orfstarts <- orfstarts[indices]
        orfstops <- orfstops[indices]
        orflengths <- numeric()
        numorfs <- length(orfstarts)
        for (i in 1:numorfs) {
                orfstart <- orfstarts[i]
                orfstop <- orfstops[i]
                orflength <- orfstop - orfstart + 1
                orflengths <- append(orflengths, orflength, after = length(orflengths))
        }
        mylist <- list(orfstarts, orfstops, orflengths)
        return(mylist)
}
max_orf <- function(OneSeq) {
        seq <- unlist(seqinr::getSequence(OneSeq,TRUE))
        orf_pos <- pos_orf(seq)
        len <- unlist(orf_pos[[3]])
        max_len <- max(len)
        if (is.na(max_len) == TRUE) {
                max_len <- 0
                max_cov <- 0
                orf_seq <- unlist(seqinr::getSequence(OneSeq, FALSE))
        } else {
                max_cov <- max_len / nchar(seq)
                max_pos <- which(len == max_len)[1]
                orf_seq <- unlist(seqinr::getSequence(OneSeq, FALSE))[orf_pos[[
                        1]][max_pos]:orf_pos[[2]][max_pos]]
        }
        res <- list(ORF = orf_seq, ORF_len = max_len, ORF_cov = max_cov)
        return(res)
}

#' Find the ORFs
#' @description This is a function to find the ORFs in one sequence.
#' @param OneSeq Is one sequence. Can be an object of the class SeqFastadna or
#' just the class character contains the sequence.
#' @return Returns a data frame. The first row is the ORF regions, the second row
#' is the lengths of the ORFs and the third row is the ORFs' coverages.
#' @author Han Siyu (Developed from the R code by Avril Coghlan)
#' @details This function can extract all ORFs of one sequence. It returns the
#' regions, lengths and coverages of ORFs. Coverage is the the ratio of ORF to
#' transcript lengths. If there are no ORF in one sequence, the first row
#' will display "No ORF is found in the sequence" while the length and coverage
#' will be zero. The package "seqinr" will be attached automatically when the
#' package "lncpred" are loaded. Users can use the function of "seqinr" to handle
#' their data.
#' This function is developed from the Avril Coghlan's R code which is used to
#' find the positions of ORFs. The original code is slightly modified here because
#' the package it depends on has updated.
#' @importFrom seqinr getSequence
#' @examples
#' ### For one sequence: ###
#' OneSeq <- c("cccatgcccagctagtaagcttagcc")
#' Seq_ORF1 <- find_orfs(OneSeq)
#' ### For a FASTA file contains several sequences: ###
#' ### Use the "read.fasta" function of package "seqinr" to read the file: ###
#' Seqs <- read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#' ### Use apply function to find ORFs: ###
#' Seq_ORF2 <- sapply(Seqs, find_orfs)
#' @export

find_orfs <- function(OneSeq) {
        seq <- unlist(seqinr::getSequence(OneSeq,TRUE))
        orf_pos <- pos_orf(seq)
        len <- unlist(orf_pos[[3]])
        if (length(orf_pos[[1]]) == 0) {
                orf_seq <- c("No ORF is found in this sequence.")
                res <- data.frame(ORF = orf_seq, ORF_length = 0,
                                  ORF_coverage = 0, stringsAsFactors = FALSE)
        } else {
                res <- data.frame()
                for (i in 1:length(len)) {
                        cov <- len[i] / nchar(seq)
                        orf_seq <- substr(seq, orf_pos[[1]][i], orf_pos[[2]][i])
                        orf <- data.frame(ORF = orf_seq, ORF_length = len[i],
                                          ORF_coverage = cov, stringsAsFactors = FALSE)
                        res <- rbind(res, orf)
                }
        }
        return(res)
}

#' Extract the Sequence Features
#' @description This is a funcion to extract the features when users want to
#' build their own model.
#' @param Seq Are the sequences that users want to extract the features. Should
#' be an object of the class SeqFastadna.
#' @param label Optional. A character which indicates the label of the sequences.
#' @param with.parallel Logical. If TRUE (Default), the process will be run in parallel.
#' @param cl.core The number of cores used to create cluster. Defualt value is all the
#' CPU cores available. (Obtain by function parallel::detectCores())
#' @return Returns a data frame contains the features. The values are numeric.
#' If uses provide a label, the values of the "label" column are factors.
#' @author Han Siyu (The extraction of ORF is based on the code by Avril Coghlan)
#' @details This is a funcion to extract the features. For each sequence, there
#' will be 1366 features which consist of the length, coverage of the longest ORF
#' and the frequencies of 1~5 adjoining-base(s) in the longest ORF region. For 1
#' base, there will be A/C/G/T, i.e. 4 features, and for 2 adjoining-bases, there
#' will be AA/AC/AG...TC/TG/TT, 4^2 features. That is to say there will be
#' 4+4^2+4^3+4^4+4^5 features of frequency. If there are more than one longest
#' ORF, only the first one will be considered. And if there is no ORF in the
#' sequence, the length and coverage of the sequence will be 0 while the
#' frequencies will be calculated on the whole sequence. Users can use the
#' data frame returned by this function to build their own model. To build a svm
#' model, users can use "svm" function of package "e1071" and please refer to its
#' documentation for further details. The package "e1071" will be loaded
#' antumatically when users employ the package "lncpred".
#' The extraction of ORF features is based on Avril Coghlan's R code. The original
#' code is to find the positions of the start and stop codons of the ORFs, and
#' the code is slightly modified because the package it depends on has updated.
#' @importFrom seqinr count
#' @importFrom e1071 svm
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @examples
#' ### Use the "read.fasta" function of package "seqinr" to read the file: ###
#' Seqs <- read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#' ### Without the label: ###
#' features_data <- extract_features(Seqs)
#' ### Label attached to every sequence: ###
#' features_data1 <- extract_features(Seqs[1:3], "label one")
#' features_data2 <- extract_features(Seqs[4:6], "label two")
#' training_set <- rbind(features_data1, features_data2)
#' ### Users can use "svm" function of package "e1071" to build a new model. ###
#' ### The label needs to be attached before training the new model. ###
#' ### This is only an example and you may see some warnings.  ###
#' new_model <- svm(label ~ ., data = training_set, scale = TRUE,
#'                  kernel = "radial", probability = TRUE)
#' test_set <- extract_features(Seqs[7])
#' pred <- predict(new_model, test_set, probability = TRUE)
#' ### For further details of function "svm" please refer to the documentation of package "e1071". ###
#' @export

extract_features <- function(Seq, label = NULL, with.parallel = TRUE, cl.core = parallel::detectCores()) {
        if(with.parallel == TRUE){
                cl <- parallel::makeCluster(cl.core)
                parallel::clusterExport(cl, varlist = c("findPotentialStartsAndStops", "pos_orf"), envir = environment())

                orfs <- parallel::parSapply(cl, Seq, max_orf)
                orf_cov <- data.frame(coverage = as.numeric(orfs[3, ]))
                orf_len <- data.frame(length = as.numeric(orfs[2, ]))

                freq1 <- data.frame(t(parallel::parSapply(cl, orfs[1, ], seqinr::count,
                                                          wordsize = 1, freq = TRUE)))
                freq2 <- data.frame(t(parallel::parSapply(cl, orfs[1, ], seqinr::count,
                                                          wordsize = 2, freq = TRUE)))
                freq3 <- data.frame(t(parallel::parSapply(cl, orfs[1, ], seqinr::count,
                                                          wordsize = 3, freq = TRUE)))
                freq4 <- data.frame(t(parallel::parSapply(cl, orfs[1, ], seqinr::count,
                                                          wordsize = 4, freq = TRUE)))
                freq5 <- data.frame(t(parallel::parSapply(cl, orfs[1, ], seqinr::count,
                                                          wordsize = 5, freq = TRUE)))
                parallel::stopCluster(cl)
        } else {
                orfs <- sapply(Seq, max_orf)
                orf_cov <- data.frame(coverage = as.numeric(orfs[2, ]))
                orf_len <- data.frame(length = as.numeric(orfs[2, ]))

                freq1 <- data.frame(t(sapply(orfs[1, ], seqinr::count, wordsize = 1, freq = TRUE)))
                freq2 <- data.frame(t(sapply(orfs[1, ], seqinr::count, wordsize = 2, freq = TRUE)))
                freq3 <- data.frame(t(sapply(orfs[1, ], seqinr::count, wordsize = 3, freq = TRUE)))
                freq4 <- data.frame(t(sapply(orfs[1, ], seqinr::count, wordsize = 4, freq = TRUE)))
                freq5 <- data.frame(t(sapply(orfs[1, ], seqinr::count, wordsize = 5, freq = TRUE)))
        }

        if(is.null(label)){
                feature_data <- cbind(orf_len, orf_cov, freq1, freq2, freq3,
                                      freq4, freq5)
        } else {
                feature_data <- cbind(label, orf_len, orf_cov, freq1, freq2,
                                      freq3, freq4, freq5)
        }
        return(feature_data)
}

#' LncRNAs Identification (Default Model)
#' @description Using the default model to identify the sequences.
#' @param Seq The sequence(s) needed to be identified. Can be an object of the
#' class SeqFastadna or the class character.
#' @param species A String indicates the species name. Use "human", "mouse" or
#' "c.elegans" to specify which model is used to predict the sequences.
#' @param detail If TRUE, the result will provide the class of the sequence as
#' well as the coding potential and the length and coverage of the longest ORF.
#' Else, only the class of the sequences will be returned.
#' @param with.parallel Logical. If TRUE (Default), the process will be run in parallel.
#' @param cl.core The number of cores used to create cluster. Defualt value is all the
#' CPU cores available. (Obtain by function parallel::detectCores())
#' @return Returns a data frame indicates the class of the sequence. The coding
#' potential and information of the longest ORF will also be provided if the
#' "detail" option is set as TRUE.
#' @author Default training model is developed by Han siyu. The model is built
#' with the LIBSVM in e1071 package.
#' @details Utilizing the default model to predict the sequences. The default
#' model for "human" is trained on human's long non-coding RNAs(lncRNAs) and
#' protein-coding transcripts, the model for "mouse" is trained on mouse's lncRNAs
#' and coding sequences(CDs); and the model for "c.elegans" is trained on non-coding
#' RNAs (ncRNAs) and CDs of Caenorhabditis elegans, this model can also be
#' applied to some invertebrata such as Saccharomyces cerevisiae.
#' @importFrom stats predict
#' @importFrom parallel detectCores
#' @examples
#' ### Use the "read.fasta" function of package "seqinr" to read the file: ###
#' Seqs <- read.fasta(file =
#' "http://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt")
#' ### Predict the Sequences: ###
#' pred1 <- lncident(Seqs, detail = TRUE)
#' ### Predict one Sequence: ###
#' pred2 <- lncident(c("cccatgcccagctagtaagcttagcc"), species = "c.elegans", with.parallel = TRUE, detail = FALSE)
#' ### Note that species name needs to be string and has no mistake. ###
#' pred3 <- lncident(c("cccatgcccagctagtaagcttagcc"), species = "C.elegans", detail = FALSE)
#' ### Will get a warning. ###
#' @export

lncident <- function(Seq, species = "human", with.parallel = TRUE, cl.core = parallel::detectCores(), detail = FALSE) {
        if(class(species) == "character" &&
           (species == "human" || species == "mouse" || species == "c.elegans")) {
                if(species == "human"){
                        model = human
                } else if(species == "mouse"){
                        model = mouse
                } else {
                        model = c.elegans
                }

                test_set <- extract_features(Seq, with.parallel = with.parallel, cl.core = cl.core)
                pred <- stats::predict(model, test_set, probability = TRUE)
                if(detail == TRUE) {
                        res <- data.frame(Prediction = pred,
                                          Coding.Potential = attr(pred, "probabilities")[ ,2],
                                          Max.ORF.Length = test_set[ ,1],
                                          Max.ORF.Coverage = test_set[ ,2])
                        return(res)
                } else {
                        res <- data.frame(Result = pred)
                        return(res)
                }
        } else {
                print("Wrong Species Name.")
        }

}
