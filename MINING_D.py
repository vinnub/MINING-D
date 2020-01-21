from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from joblib import Parallel, delayed
import networkx as nx
import click

import numpy as np
from scipy.stats import chisquare

########################################################################################################################
########################################  Some Utility functions ######################################################
########################################################################################################################

def get_extension_stats(initial_mer, direction, cdr3s):
    '''
    Get counts for the 4 possible extensions corresponding to either left or right extensions of 'initial_mer' and the chi-square p-value. 
    Returns (p_value, extended_mers, counts_of_extended). 
    '''
    assert direction in ['left', 'right'], 'Direction must be either left or right'
    if direction == 'right':
        extended_mers = [initial_mer + "A", initial_mer + "G",initial_mer + "C",initial_mer + "T" ]
    else: 
        extended_mers = ["A" + initial_mer ,"G" +  initial_mer ,"C" + initial_mer ,"T" + initial_mer  ]
    
    count_of_extended = _abundances_of(extended_mers, cdr3s) 
    p_value = chisquare(list(count_of_extended.values()))[1] if sum(count_of_extended.values())>0 else 1
    return p_value, count_of_extended

def _abundances_of(sequences, seq_list):
    ''' checks the abundance of 4 extended mers '''
    abundances  = {}
    for seq_to_check in sequences:
        abundances [seq_to_check] = 0
        for a_seq in seq_list: 
            abundances [seq_to_check] += int( seq_to_check in a_seq)
    return abundances 

def extensions_for_abundant_kmers(a_CDR3_rep_inst, k_mer_len, num_k_mers_to_extend, alpha, n_cores):
    '''
    Returns a dictionary where keys are the 'num_k_mers_to_extend' most abundant k-mers in the 'a_CDR3_rep_instance'. 
    The value for a key are all the extensions corresponding to the key. The computation is done parallely using n_cores. 
    '''
    seqs_to_extend = a_CDR3_rep_inst.get_abundant_k_mers(k_mer_len = k_mer_len, num=num_k_mers_to_extend, )
    extensions_list = Parallel(n_jobs = n_cores)(delayed(a_CDR3_rep_inst.extensions_of)(a_seq, alpha= alpha, k_mer_len=k_mer_len) for a_seq in seqs_to_extend)
    return dict(zip(seqs_to_extend, extensions_list))

########################################################################################################################
######################################## CDR3 repertoire class #########################################################
########################################################################################################################

class CDR3_rep:
    '''
    Variables needed= k, num_k_mers, p_value_threshold
    
    '''
    
    def __init__ (self, filename):
        def get_unique_cdr3s(a_cdr_file):
            return {str(records.seq) for records in SeqIO.parse(a_cdr_file, "fasta")}
        
        self.filename = filename
        self.cdr3s = get_unique_cdr3s(self.filename)
    
    # Find the mean relative position in CDR database 
    def rel_pos_in_cdr3s(self, seq):
        rel_pos = list()
        for cdr in self.cdr3s:
            if seq in cdr: rel_pos.append(cdr.index(seq)/(len(cdr)-len(seq)+1))
        return np.mean(rel_pos)

    def get_abundant_k_mers(self, k_mer_len , num ):
        '''
        Return 'num' most abundant k-mers
        k_mer_length : Length of k-mer
        '''
        abundances = {}
        for a_cdr3 in self.cdr3s:
            for i in range(len(a_cdr3) - k_mer_len + 1):
                k_mer = a_cdr3[i:i+k_mer_len]
                if k_mer not in abundances:
                    abundances[k_mer] =1
                else:
                    abundances[k_mer] += 1
        sorted_abundances_list = sorted(abundances.items(), key=lambda item: item[1], reverse = True)
        abundant_k_mers = [k_mer for k_mer, abundance in sorted_abundances_list[:num]]
        return abundant_k_mers

    def _get_extnsn_and_suggstns(self, initial_mer, alpha , verbose = False, ):
        ''' Get one extension from the initial_mer and other sequences to extend. 
        Returns (one_extension (a string) , other_sequences_to_extend (list of strings))'''
        initial_mer_copy = initial_mer
        p_value_right = 0
        p_value_left = 0
        p_value_3 = 1
        suggested_extensions = list()
        initial_mer_cdr3s = {a_cdr3 for a_cdr3 in self.cdr3s if initial_mer in a_cdr3}

        while (min(p_value_right, p_value_left) < alpha):
            p_value_right, count_of_extended_right = get_extension_stats(initial_mer, 'right', initial_mer_cdr3s)
            p_value_left, count_of_extended_left  = get_extension_stats(initial_mer, 'left', initial_mer_cdr3s)

            #Print the counts on the left and right 
            if verbose:
                print(f"For Right extension {count_of_extended_right} \nThe p-value is {p_value_right}\n" )
                print(f"For Left extension {count_of_extended_left} \nThe p-value is {p_value_left}\n" )

            if min(p_value_right, p_value_left) < alpha:
                count_of_extended = count_of_extended_right if p_value_right <= p_value_left else count_of_extended_left
                all_max = [key for key, value in count_of_extended.items() if value == max(count_of_extended.values())]

                if len(all_max) == 1: # if there is only 1 max in counts list
                    small_counts = {k:v for k,v in count_of_extended.items() if k != all_max[0]}
                    if sum(small_counts.values()) >0:
                        p_value_3 = chisquare(list(small_counts.values()))[1]
                        if verbose: print(f"The remaining counts are: {small_counts} \nThe p-value is {p_value_3}\n" )

                        if p_value_3 < alpha:
                            suggested_extensions.append(sorted(small_counts.items(), key = lambda item : item[1], reverse = True)[0][0])
                            if verbose: print("\nAlso look for extension of %s" %suggested_extensions[-1] )        

                        initial_mer = all_max[0]
                        if verbose: print("Add %s to the right \n\n" %initial_mer[-1]) if p_value_right <= p_value_left else print("Add %s to the left \n\n" %initial_mer[0])
                    else:
                        initial_mer = all_max[0]
                        if verbose: print("Add %s to the right \n\n" %initial_mer[-1]) if p_value_right <= p_value_left else print("Add %s to the left \n\n" %initial_mer[0])

                else: # if there are multiple max in counts list
                    initial_mer = all_max[0]
                    if verbose: print("Add %s to the right \n\n" %initial_mer[-1]) if p_value_right <= p_value_left else print("Add %s to the left \n\n" %initial_mer[0])
                    for i in range(1,len(all_max)):
                        suggested_extensions.append(all_max[i])
                        if verbose: print("\nAlso look for extension of %s" %suggested_extensions[-1] )        

        return((initial_mer, suggested_extensions))

    def get_multiple_extensions(self, initial_mer, alpha, verbose = False):
        '''
        Returns multiple extensions starting from one sequence.'''
        extensions = list()
        confirmed_extension, suggested_extensions = self._get_extnsn_and_suggstns(initial_mer, alpha = alpha, verbose = verbose)
        extensions.append(confirmed_extension)
        if len(suggested_extensions)>0:
            index = 0
            while(index < len(suggested_extensions)):
                a_suggestion = suggested_extensions[index]
                confirmed_extension, further_suggestions = self._get_extnsn_and_suggstns(a_suggestion, alpha = alpha,verbose = verbose)
                extensions.append(confirmed_extension)
                if len(further_suggestions)>0: suggested_extensions = suggested_extensions + further_suggestions
                index += 1
        return extensions

    #Shrink sequences 
    def _shrink_seq(self,  starting_seq, alpha , k_mer_len, verbose = False):
        seq = starting_seq
        while 1:
            if len(seq) < k_mer_len : break

            truncated_seq_left = seq[1:]
            truncated_seq_right = seq[:-1]

            p_value_left, count_of_extended_left  = get_extension_stats(truncated_seq_left, 'left', self.cdr3s)
            p_value_right, count_of_extended_right = get_extension_stats(truncated_seq_right, 'right', self.cdr3s)

            if sorted(count_of_extended_left.items(), key = lambda item : item[1], reverse = True)[0][0]!= seq:
                if verbose ==True: print ("Delete leftmost letter")
                seq = truncated_seq_left
                if verbose ==True: print("The remaining sequence is %s" %seq)
                continue

            if sorted(count_of_extended_right.items(), key = lambda item : item[1], reverse = True)[0][0]!= seq:
                if verbose ==True: print ("Remove right-most letter")
                seq = truncated_seq_right
                if verbose ==True: print("The remaining sequence is %s" %seq)
                continue

            if max(p_value_right, p_value_left) > alpha:
                if verbose ==True: print("Delete right letter") if p_value_right >= p_value_left else print ("Delete left letter ") 
                seq =   truncated_seq_right if p_value_right >= p_value_left else  truncated_seq_left
                if verbose ==True: print("The remaining sequence is %s" %seq)
            else:
                break 
        return seq

    def extensions_of(self, a_seq, alpha, k_mer_len ):
        '''
        Get all processed extensions of the starting sequence 'a_seq'. 
        Processed here refers to shrinking the extensions and then removing extensions that are subsequences of other extensions.
        ''' 
        def shrink_and_remove_subseqs(ext_list, alpha):
            a_list = list(np.unique(ext_list))
            shrunk_list = [self._shrink_seq(item, alpha = alpha, k_mer_len = k_mer_len) for item in a_list]
            shrunk_lens = [len(item) for item in shrunk_list]
            sorted_list = [shrunk_list[i] for i in np.argsort(shrunk_lens)]
            to_remove = list()
            for i in range(len(sorted_list)):
                for j in range((i+1), len(sorted_list)):
                    if sorted_list[i] in sorted_list[j]:
                        to_remove.append(sorted_list[i])
                        break
            final_list = list (set(sorted_list).difference(set(to_remove)))  
            return final_list

        unprocessed_extensions = self.get_multiple_extensions(a_seq, alpha=alpha)
        return shrink_and_remove_subseqs(unprocessed_extensions, alpha = alpha)
    
########################################################################################################################
#################################### Processing extensions and candidate Ds   ##########################################
########################################################################################################################

def filter_unidirectional_extensions(an_extensions_dict, side_alpha):
    '''
    Takes as input a dictionary of extensions where keys are kmers and values are extensions and filters 
    those extensions that are unidirectional according to the definition given in the paper.
    Returns a new dictionary and does not change the original. 
    '''
    def keep_bidirectional_extensions (starting_seq, extensions):
        '''
        Takes as input a kmer and its extensions and returns only extensions that are unidirectional w.r.t. the kmer'''
        bidir_ext = set()
        for ext in extensions:
            alignment = pairwise2.align.localms(starting_seq, ext , 1, -100, -100, 0, penalize_end_gaps = False)[0]
            left_ext = alignment[3]
            right_ext = len(alignment[0]) - alignment[4]
            if alignment[0][0] == '-' and alignment[0][-1] == '-':
                if ((abs(left_ext-right_ext))/max(left_ext, right_ext))  < side_alpha:
                    bidir_ext.add(ext)
        return bidir_ext

    filtered = {}
    for kmer, extensions in an_extensions_dict.items():
        filtered[kmer] = keep_bidirectional_extensions(kmer, extensions)
    return filtered


def graph_from_list(a_list, threshold = 1):
    '''
    This function returns a networkx graph from a list, where there is an edge between two nodes if and only if the 
    score between two corresponding extensions in the list is within the specified threshold
    '''
    score_matrix = np.zeros((len(a_list), len(a_list)))
    for i in range(len(a_list)):
        for j in range(i + 1,len(a_list)):
            alignment = get_alignment(a_list[i], a_list[j] )
            score_matrix[i,j] = max(min(len(a_list[i]), len(a_list[j])),8) - alignment[2]
    
    score_matrix[score_matrix > threshold] = 0
    score_matrix[score_matrix > 0] = 1

    a_graph = nx.from_numpy_matrix(score_matrix)
    return a_graph

######################################## Merging Cliques ######################################################

def merge_cliques(a_list, a_cdr3_rep, alpha, k_mer_len, graph_threshold =2, verbose = False):
    '''
    This function merges the cliques in a graph_from_list formed from a_list with the paramter graph threshold. They are 
    merged in the following way:
    1. Take the intersection of the largest clique in all connected components. 
    2. Extend the intersection. Make sure the extensions are central in CDRs. 
        In this version of the function, step3 is not performed. We include all the extensions which we get from an intersection. 
        3. If the extensions from one intersection are within a score of extensions_threshold, run remove_score_s on those.
        Step 3 is included in merge_cliques
    4. Repeat for all connected components. 
    '''
    
    def _align_intersection(indices, a_list, verbose = False):
        ''' 
        Function that takes into input a list of indices and a list of sequences, and returns the intersection of the 
        pairwise alignments 
        ''' 
        middle_extensions = a_list
        alignment = get_alignment(middle_extensions[indices[0]],middle_extensions[indices[1]] )  
        intersection = alignment[0][alignment[3]:alignment[4]]
        if verbose ==True: print( "The first intersection is %s" %intersection)
        for i in range(2, len(indices)):
            alignment = get_alignment(intersection,middle_extensions[indices[i]] )   
            intersection = alignment[0][alignment[3]:alignment[4]]
            if verbose ==True: print("The %s th intersection is %s" %(i,intersection) )
        return intersection 

    to_remove = list()
    to_add= list()
    graph = graph_from_list(a_list, graph_threshold )
    subgraphs  = list(nx.connected_component_subgraphs(graph))
    for i, a_subg in enumerate(subgraphs):
            if a_subg.number_of_nodes() == 1: 
                if verbose: print("Only 1 node in subgraph %d" %i)
            else:
                largest_clique_nodes =  sorted(list(nx.find_cliques(a_subg)), key = len)[-1]
                if verbose: print(largest_clique_nodes)
                intersection = _align_intersection(largest_clique_nodes, a_list)
                to_remove.extend([a_list[index] for index in largest_clique_nodes ])
                extensions = a_cdr3_rep.extensions_of(intersection, alpha = alpha, k_mer_len = k_mer_len)
                ext_in_middle = [ext for ext in extensions if (a_cdr3_rep.rel_pos_in_cdr3s(ext)> 0.2) and (a_cdr3_rep.rel_pos_in_cdr3s(ext)<0.7) ]
                to_add.extend(ext_in_middle)
                if verbose: print("Removing %d  Adding %d" %(len(to_remove), len(to_add)))
               
    return list ( (set(a_list).difference(set(to_remove))  ).union(set(to_add)))

def merge_cliques_iteratively(a_list, a_cdr3_rep, alpha, k_mer_len, graph_threshold = 2,  verbose = False):
    '''
    This function runs the function merge_cliques iteratively until the lengths of the list in two consecutive 
    iterations are same. Returns a new list. 
    '''
    most_recent = a_list
    initial_length = len(most_recent)
    i = 1
    while(1):
        if verbose: print("Doing %d iteration" %i)
        temp_list = (merge_cliques(most_recent, a_cdr3_rep=a_cdr3_rep, alpha=alpha, k_mer_len=k_mer_len,
                                   graph_threshold =graph_threshold, verbose = verbose))
        #delete extra Gs
        temp_db = D_gene_Candidates()
        temp_db.add_sequences(temp_list)
        temp_db_new = temp_db.del_abund_rand_insrtns()
        most_recent = temp_db_new.sequences #a list
        if verbose: print("Sequences went from %d to %d" %(initial_length, len(most_recent)))
        if len(most_recent) == initial_length:
            break
        else:
            initial_length = len(most_recent)
            i += 1
    return most_recent

################################# Final Processing (Merge similar candidates etc.) #######################################

def _action (seq1, seq2, verbose = False):
    '''
    This function tells the operation to be applied to the two sequences that are deemed to be 
    close to one another in the previous step. 
    (usually if there is an edge from the function graph_from_list)

    If - is only in 1 sequence, keep the bigger sequence.
    If - is not present in any sequence, keep the intersection. 
    If - is present in both, see if the sequences can be merged. If there are nucleotides that 
    contradict in the two sequences, do not merge. If there are no contradicting nucleotides, merge.
    ''' 
    align = pairwise2.align.localms(seq1, seq2 , 1, -100, -100, 0, penalize_end_gaps = False)[0]
    if verbose == True: print(align)
    first_seq, second_seq, score, start, end = align
    if verbose == True: print(format_alignment(*align))
   
    if ('-' in first_seq) and ('-' in second_seq):
        to_check = [item for item in range(len(first_seq)) if item not in range(start, end)]
        same_nuc = [(first_seq[letter] != '-') and (second_seq[letter] != '-') for letter in to_check]
        
        if any(same_nuc):        
            #if verbose == True: print('Unmatching letter at position %d' %letter) 
            if verbose == True: print('Keep the bigger sequence')
            return seq1 if len(seq1)>= len(seq2) else seq2
        else:
            if verbose == True: print('Merge sequences')
            return _merge_sequences(first_seq, second_seq, score, start, end)
    
    elif ('-' in first_seq) or ('-' in second_seq):
        if verbose == True: print('Keep the bigger sequence')
        return seq1 if len(seq1)>= len(seq2) else seq2
    
    else:
        if verbose == True: print('Keep the intersection')
        return(seq1[start:end])

    

def _merge_sequences(first_seq, second_seq, score, start, end):
    '''
    This function merges two sequences. If one sequence has '-' in the alignment and the other 
    has some nucletotide, it keeps the nucleotide. 
    --ACCTCTT
    TCACCTC-- will return 

    TCACCTCTT.
    This is only done for seuqnces that were deemed appropriate to merge in the function '_action'.
    '''
    part_one = first_seq[:start] if all([letter == '-' for letter in second_seq[:start] ]) else second_seq[:start]
    part_two = first_seq[start:end]
    part_three = first_seq[end:] if all([letter == '-' for letter in second_seq[end:] ]) else second_seq[end:]
    
    return ''.join([part_one, part_two, part_three])
    

    
def process_edges_in_graph(d_list, threshold = 2, verbose = False):
    '''
    Draw a graph from a list. Start taking '_actions' on the nodes, edge by edge. 
    If any edges is processed, skip all the other edges containing the nodes in the already processed
    edges. 
    This can also be done iteratively if needed. 
    '''
    g1 = graph_from_list(d_list, threshold=threshold)
    edges = list(g1.edges)
    done_nodes = []
    to_remove = []
    to_add = []
    for edge in edges:
        if verbose ==True: 
            print('\n')
            print("DONE NODEs", done_nodes)
        if any([node in done_nodes for node in edge]): 
            if verbose ==True: print('Skipping edge ' ,edge)
            continue
        else:
            done_nodes.extend(edge)
            if verbose ==True: print('Processing for edge' , edge)
            to_remove.extend([d_list[node] for node in edge])
            add = _action((*[d_list[node] for node in edge]), verbose=verbose)
            if verbose == True: print(add)
            to_add.append(add)
    return list ( (set(d_list).difference(set(to_remove))  ).union(set(to_add)))



########################################################################################################################
######################################## D genes database class ######################################################
########################################################################################################################

##### Utility functions. 

def get_alignment(seq1, seq2):
    return pairwise2.align.localms(seq1, seq2 , 1, -100, -100, 0, penalize_end_gaps = False)[0]  

def score(seq1, seq2, verbose = False):
    '''Function that calculates the score metric between two sequences '''
    alignment = get_alignment(seq1, seq2)
    if verbose ==True: 
        print(format_alignment(*alignment))
    return max(min(len(seq1), len(seq2)) ,8) - alignment[2]

def extra_symbols_in_smaller(alignment, check_seq2 = False):
    '''Takes as input an alignment and gives the extra symbols present in the seq1 (or seq2 if check_seq2 == True).'''
    if not check_seq2:
        return {alignment[0][i] for i in set(range(len(alignment[0]))).difference(set(range(alignment[3],alignment[4])))}
    else:
        return {alignment[1][i] for i in set(range(len(alignment[0]))).difference(set(range(alignment[3],alignment[4])))}
    
    
##### Class definitions. 
class d_segment:
    '''Represents a D gene segment.'''
    def __init__(self, name = '', seq = ''):
        self.name = name
        self.seq = seq
    def __str__(self):
        to_print = ''
        for k,v  in self.__dict__.items():
                to_print += str(k) + ' ' + str(v) + '\n'
        return to_print
    

class D_gene_Candidates:
    '''A list of candidate d segments'''
    def __init__(self, name = '', segments = None ): 
        #segments is a list of d_segments
        self.name = name
        self.segments = segments if segments is not None else []
        
    def __str__(self):
        to_print = ''
        for a_seg in self.segments:
            for k,v  in a_seg.__dict__.items():
                to_print += str(k) + ' ' + str(v) + '\n'
            to_print += '\n'
        return to_print
    
    def write_to_fasta(self, file_name):
        new_list = list()
        for seg in self.segments:
            new_list.append( SeqRecord(Seq(seg.seq),id= '%s'%seg.name,description="" ))
        SeqIO.write(new_list, file_name, "fasta")

    @property
    def sequences(self):
        return [str(d_seg.seq) for d_seg in self.segments]
    
    @property
    def num_seq(self):
        return len(self.segments)
    
    @property
    def seg_names(self):
        return [str(d_seg.name) for d_seg in self.segments]
    
    def add_d_segment(self, d_seg):
        self.segments.append(d_seg)
    
    def add_sequences(self, a_seq_list):
        for i, a_seq in enumerate(a_seq_list):
            self.add_d_segment(d_segment(seq= a_seq, name = 'Cand_' + str(i+1)))
    
    def add_RP_in_cdrs_frm(self, cdr3_rep):
        for a_dseg in self.segments:
            a_dseg.RP = cdr3_rep.rel_pos_in_cdr3s(a_dseg.seq)
            
    def get_central_d_segments(self, left_th = 0.2, right_th = 0.7):
        central_d_segs = []
        for a_dseg in self.segments:
            if (a_dseg.RP > left_th) and (a_dseg.RP < right_th):
                central_d_segs.append(a_dseg)
        return D_gene_Candidates(segments=central_d_segs, name = 'central')
    
    def filter_cand_with_motif(self, motif = 'TACTACTAC'):
        new_list = [a_seg.seq for a_seg in self.segments if motif not in a_seg.seq]
        to_return_db = D_gene_Candidates()
        to_return_db.add_sequences(new_list)
        return to_return_db

    def merge_cliques_in_db(self, a_cdr3_rep, alpha, k_mer_len, graph_threshold = 2,  verbose = False):
        '''
        This method merges the cliques iteratively in the graph on the sequences and returns a new D_gene_Candidates object.
        '''
        new_seqs = merge_cliques_iteratively(self.sequences, a_cdr3_rep = a_cdr3_rep, 
                                             alpha = alpha, k_mer_len = k_mer_len,
                                             graph_threshold=graph_threshold, 
                                             verbose = verbose)
        to_return_db = D_gene_Candidates()
        to_return_db.add_sequences(new_seqs)
        return to_return_db

    def merge_smlr_seqs(self, graph_threshold = 2,  verbose = False):
        '''
        This method does the final processing on the candidates (merge very similar sequences etc.) returns a new D_gene_Candidates object.
        '''
        new_seqs = process_edges_in_graph(self.sequences, threshold=graph_threshold, verbose = verbose)
        to_return_db = D_gene_Candidates()
        to_return_db.add_sequences(new_seqs)
        return to_return_db
    
    def del_abund_rand_insrtns(self, symbol_to_del= 'G', g_threshold = 3, verbose = False): 
        ''' 
        This function returns a new D_gene_candidates object with some segments deleted from the original object. 
        A segment that is at least a score of g_threshold from another segment of atleast the same length, and the score is only due 
        to the 'symbol_to_del', is deleted. 
        eg. AAAAAAAACTCTC
            AAAAAAAAGGG
        In the above example the second one would be deleted, because it's been observed that random insertions of G are 
        very common. 
        '''

        a_list = self.segments
        lengths = [len(item.seq) for item in a_list]
        sorted_list = [a_list[i] for i in np.argsort(lengths)]
        to_remove = list()
        for i in range(len(sorted_list)):
            for j in range(i + 1, len(sorted_list)):
                alignment = get_alignment(sorted_list[i].seq, sorted_list[j].seq)
                s = score(sorted_list[i].seq, sorted_list[j].seq)
                if s >= g_threshold: 
                    if len(sorted_list[j].seq) > len(sorted_list[i].seq):
                        if len(extra_symbols_in_smaller(alignment).difference({'-',symbol_to_del})) == 0:
                            to_remove.append(sorted_list[i])
                            if verbose: 
                                print("Removing %s because of %s" %(sorted_list[i].seq, sorted_list[j].seq))
                                score(sorted_list[i].seq, sorted_list[j].seq, verbose = True)
                            break
                    else:
                        alignment0 = len(extra_symbols_in_smaller(alignment).difference({'-',symbol_to_del}))
                        alignment1 = len(extra_symbols_in_smaller(alignment, check_seq2=True).difference({'-',symbol_to_del}))
                        if alignment0 == 0 or alignment1 ==0:
                            to_remove.append(sorted_list[i]) if alignment0 ==0 else to_remove.append(sorted_list[j])
                            if verbose: 
                                print("Removing %s because of %s" %(sorted_list[i].seq, sorted_list[j].seq)) if alignment0 ==0 else print("Removing %s because of %s" %(sorted_list[j].seq, sorted_list[i].seq))
                                score(sorted_list[i].seq, sorted_list[j].seq, verbose = True)
        new_list = list (set(sorted_list).difference(set(to_remove))  )
        to_return = D_gene_Candidates(segments=new_list)
        return to_return

    
    
@click.command() 
@click.option('-i', '--input_cdr_file',  required=True, help='Input file with consensus CDR3s')
@click.option('-o', '--output_file', required=True, help='Output file to store inferred D genes')
@click.option('-k',  default = 10,  show_default=True, help='Length of starting k-mers')
@click.option('-n', '--num_k_mers', 'k_mers',  default = 300, show_default=True, help='Number of most abundant k-mers to extend')
@click.option('-p','--p_val_th',  default = 4.5*10**(-36), show_default=True, help='p-value threshold for the extension procedure')
@click.option('-t', '--n_cores', default = 10, show_default=True,help='Number of cores')
@click.option('--bidir_alpha',  default = 0.5, show_default=True,help='Alpha for filtering bidirectional extensions (see paper)')
@click.option('-g', '--cliq_th', default = 2, show_default=True,help='Similarity metric threshold for generating graphs')
def _run_MINING_D(input_cdr_file,output_file,  k,k_mers, p_val_th, n_cores, bidir_alpha, cliq_th):
    
    
    cdr3s = CDR3_rep(input_cdr_file)
    
    # Compute multiple extensions
    extensions_dict = extensions_for_abundant_kmers(cdr3s, 
                                                    num_k_mers_to_extend= k_mers, 
                                                    n_cores=n_cores, 
                                                    k_mer_len=k, 
                                                    alpha= p_val_th)
    filtered_extensions = filter_unidirectional_extensions(extensions_dict, side_alpha=bidir_alpha)
    
    # Get a candidate set
    candidate_seqs = {ext for kmer in filtered_extensions for ext in filtered_extensions[kmer]}

    # Create a D_gene_Candidates object 
    cand = D_gene_Candidates()
    cand.add_sequences(candidate_seqs)
    
    # Add relative positions and keep only central candidates
    cand.add_RP_in_cdrs_frm(cdr3s)
    cand_central = cand.get_central_d_segments() 
    
    # Filtering using merging cliques etc. 
    cand_central = cand_central.del_abund_rand_insrtns()
    cand_after_merging_cliques = cand_central.merge_cliques_in_db(a_cdr3_rep=cdr3s, k_mer_len=k, alpha= p_val_th,
                                                                  graph_threshold=cliq_th)
    cand_after_merging_cliques = cand_after_merging_cliques.\
                                del_abund_rand_insrtns(symbol_to_del = 'G',g_threshold=2).\
                                del_abund_rand_insrtns(symbol_to_del = 'C', g_threshold=2)

    # Additional filtering for special cases (rare)
    final_candidates = cand_after_merging_cliques.merge_smlr_seqs()
    final_candidates = final_candidates.filter_cand_with_motif(motif = 'TACTACTAC')
    
    # Write to Fasta 
    final_candidates.write_to_fasta(output_file)
    return final_candidates.num_seq

if __name__ == '__main__':
    _run_MINING_D()