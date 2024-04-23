import sys
amino_acid_values = {
    "E": [1.53, 0.26],  #Glu
    "A": [1.45, 0.97],  #Ala
    "L": [1.34, 1.22],  #Leu
    "H": [1.24, 0.71],  #His
    "M": [1.20, 1.67],  #Met
    "Q": [1.17, 1.23],  #Gln
    "W": [1.14, 1.19],  #Trp
    "V": [1.14, 1.65],  #Val
    "F": [1.12, 1.28],  #Phe
    "K": [1.07, 0.74],  #Lys
    "I": [1.00, 1.60],  #Ile
    "D": [0.98, 0.80],  #Asp
    "T": [0.82, 1.20],  #Thr
    "S": [0.79, 0.72],  #Ser
    "R": [0.79, 0.90],  #Arg
    "C": [0.77, 1.30],  #Cys
    "N": [0.73, 0.65],  #Asn
    "Y": [0.61, 1.29],  #Tyr
    "P": [0.59, 0.62],  #Pro
    "G": [0.53, 0.81]   #Gly
    }
sequence_given = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLEPSQDL"

List_Containing_Possible_Helices = ["-" for i in range(len(sequence_given))]
List_Containing_Possible_Betas = ["-" for i in range(len(sequence_given))]
output = ["-" for i in range(len(sequence_given))]


# Returns score of the string entered
def Helix_Score_using_P (residue) :    
    score = 0
    for aminoAcid in residue:
        for key in amino_acid_values :
            if (key == aminoAcid):
                score+= amino_acid_values[aminoAcid][0]
                # 0th index contains P(h)
    return score

# Returns score of the string entered for beta
def Beta_Score_using_P (residue) :    
    score = 0
    for aminoAcid in residue:
        for key in amino_acid_values :
            if (key == aminoAcid):
                score+= amino_acid_values[aminoAcid][1]
                # 1st index contains P(s)
    return score

def calculate(i, type) : 
    # Pass "S"/"H" on the basis of the analysis been 
    # done, i.e. S in case of Strand/Beta analysis and 
    # H in case of Helix analysis.
    # Calculates whether given residue is a valid nucleation site.
    count = 0
    # By default, assuming the case of helix
    value_to_be_subtracted = 0
    index = 0
    if (type == "S"):
        value_to_be_subtracted = 1
        index = 1
    
    for k in range (6 - value_to_be_subtracted):
        if (amino_acid_values[sequence_given[i + k]][index] >= 1) :
            count += 1
    return count

def extend_left(j, type):
    index = 0
    if (type == "S"):
        # Case of strand
        index = 1
    while j <= len(sequence_given) - 4:
        count = 0
        for k in range(0, 4):
            count = count + amino_acid_values[sequence_given[j + k]][index]
        if count >= 4:
            if (type == "S") :
                for k in range(0, 4):
                    List_Containing_Possible_Betas[j + k] = type
                j = j + 1
            elif (type == "H") : 
                for k in range(0, 4):
                    List_Containing_Possible_Helices[j + k] = type
                j = j + 1
        else:
            break

def extend_right (j, type) :
    # By default, assuming the case of helix
    index = 0
    if (type == "S"):
        # Case of strand
        index = 1
    while j - 3 >= 0:
        count = 0
        for k in range(0, 4):
            count = count + amino_acid_values[sequence_given[j - k]][index]
        if count >= 4:
            if (type == "S") :
                for k in range(0, 4):
                    List_Containing_Possible_Betas[j - k] = type
                j = j - 1
            elif (type == "H") :
                for k in range(0, 4):
                    List_Containing_Possible_Helices[j - k] = type
                j = j - 1
        else :
            break

# --------------------------------------------
# Processing the sequence recursively, and respectively extending left/right.
def process_sequence(sequence, target, length):
    i = 0
    while i < len(sequence) - length + 1:
        count = calculate(i, target)
        if count >= length - 2:
            for k in range(length):
                sequence[i + k] = target
            extend_right(i + length - 3, target)
            extend_left(i + 2, target)
        i += 1

process_sequence(List_Containing_Possible_Helices, "H", 6)
process_sequence(List_Containing_Possible_Betas, "S", 5) 
# ----------------------------------------------


# Resolving Conflicts :
# Making a dictionary, which contains key as the index, and value as the list of charachters on that index in the 2 lists

dict_conflicts_resolve = {}
Conflicting_Regions = []
resolved_List = [0]*len(sequence_given)

for i in range(len(sequence_given)):
    dict_conflicts_resolve[i] = [List_Containing_Possible_Helices[i], List_Containing_Possible_Betas[i]]
# print(dict_conflicts_resolve)


# Updating the resolved_List on the basis conflicts resolves.
i = 0
while (i<len(sequence_given)): 
    # print (i, end = ", ")
# for i in dict_conflicts_resolve:
    # Case of no conflict :
    [first, second] = dict_conflicts_resolve[i]

    if (first == "-" and second == "-") :
        # Case when neither helix nor beta
        resolved_List[i] = "-"
        i+=1

    elif (first == "-" and second == "S") :
        # Case when Beta dominates
        resolved_List[i] = "S"
        i+=1
    
    elif (first == "H" and second == "-") :
        # Case when Helix dominates
        resolved_List[i] = "H"
        i+=1

    elif (first == "H" and second == "S") : 
        array1 = []
        array2 = []
        
        pointer1, pointer2 = i, i

        while (dict_conflicts_resolve[pointer1][0] == "H" and dict_conflicts_resolve[pointer1][1] == "S"):
            array1.append(sequence_given[pointer1])
            array2.append(sequence_given[pointer2]) 
            pointer1 += 1
            pointer2 += 1
            if (pointer1 == len(sequence_given)):
                break
        Conflicting_Regions.append((i, pointer1))
        seq1 = "".join(array1)
        seq2 = "".join(array2)
        # Finding the scores by joining array1 and array2 respectively; and then comparing which 
        # charachter is dominant.

        score_helix = Helix_Score_using_P(seq1)
        score_beta = Beta_Score_using_P(seq2)
        
        if (score_helix > score_beta) : 
            for index in range(i, pointer1) : 
                resolved_List[index] = "H"
        elif (score_helix < score_beta) : 
            for index in range(i, pointer1) : 
                resolved_List[index] = "S"
        else : 
            # If the score of both the common sequences is coming out to be same, then 
            # by default giving priority to Beta over Helix.
            for index in range(i, pointer1) : 
                resolved_List[index] = "S"

        # Modify the value of i, so as to by pass the already upadated sequences.
        i = pointer1  # Can be pointer2 as well.

    else : 
        error_message = "An error has occurred while resolving conflicts!"
        sys.stdout.write(error_message)
    # print(i)

print("Considering cases only of Helical Regions : \n", "".join(List_Containing_Possible_Helices), end = "\n\n")
print("Considering cases only of Beta Regions : \n", "".join(List_Containing_Possible_Betas), end = "\n\n")

print("Portions in conflict : \n", Conflicting_Regions, end = "\n\n")
print("Final Answer after Resolving Conflicts : \n", "".join(resolved_List), end = "\n\n")
