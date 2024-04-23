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

def return_possibilities (sequence) : 
    dict_with_residue_and_indices = {}
    possibilities = []
    for i in range(len(sequence)-5):
        possibilities.append(sequence[i : 6+i])
        dict_with_residue_and_indices[i] = [i, 6+i -1]         # i represents the nth index in the possibilities list.
    return possibilities, dict_with_residue_and_indices

# Returns score of the string entered
def Helix_Score_using_P (residue) :    
    score = 0
    for aminoAcid in residue:
        for key in amino_acid_values :
            if (key == aminoAcid):
                score+= amino_acid_values[aminoAcid][0]
                # 0th index contains P(h)
    return score

def Check_Valid_Helix (residue):
    count = 0
    for aminoAcid in residue:
        for key in amino_acid_values :
            if (key == aminoAcid):
                if(amino_acid_values[aminoAcid][0] >= 1) :
                    count +=1
    if (count >= 4) : 
        return True
    return False

# sequence_given = "WHGCITVYWMTV"
sequence_given = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLEPSQDL"

residueList, dict_with_residuePosn_and_indices = return_possibilities(sequence_given)

List_Of_Helical_Regions = []
List_Containing_Possible_Helices = ["-" for i in range(len(sequence_given))]


for position in dict_with_residuePosn_and_indices:
    start, end = position, dict_with_residuePosn_and_indices[position][1] # Both included.
    # Case when first residue => neglecting the effect of front pointer.
    if (position == 0): 
        # print ("Check 1")
        # Primary check
        if (Check_Valid_Helix(residueList[position])): 
            # print ("Check Valid 1")
            # Converting "-" to "H"
            for j in range (start, end+1):
                List_Containing_Possible_Helices[j] = "H"

            i=5+1
            current_residue = residueList[position][-3:] + sequence_given[i]
            # Assuming residue to be helical, hence assigning.
            Helical_Region = residueList[position]
            while (i<len(sequence_given)):
                if (Helix_Score_using_P(current_residue) >= 4):
                    # print ("Check 2")
                    Helical_Region += current_residue[-1]
                    if (i+1 != len(sequence_given)):
                        List_Containing_Possible_Helices[i+1] = "H"
                        current_residue = current_residue[-3:] + sequence_given[i+1]
                    i+=1
                else:         
                    break

            List_Of_Helical_Regions.append(Helical_Region)

 # Case last residue => neglecting the effect of end pointer.
    elif (position == len(residueList)-1) :
        # print ("Check 4")
        # Primary check
        if (Check_Valid_Helix(residueList[position])): 
            # print ("Check Valid 2")
            # Converting "-" to "H"
            for j in range (start, end+1):
                List_Containing_Possible_Helices[j] = "H"

            i = dict_with_residuePosn_and_indices[len(residueList)-1][0] - 1
            current_residue =  sequence_given[i] + residueList[position][ : 3]
            # Assuming residue to be helical, hence assigning.
            Helical_Region = residueList[position]
            while (i>=0):
                if (Helix_Score_using_P(current_residue) >= 4):
                    # print ("Check 5")
                    Helical_Region = current_residue[0] + Helical_Region
                    List_Containing_Possible_Helices[i] = "H"
                    i-=1
                    current_residue =  sequence_given[i] + current_residue[ : 3]
                else : 
                    break
            List_Of_Helical_Regions.append(Helical_Region)

    # Case of residues other than 1st and last.
    else :
        # print ("Check 7")
        # Primary check
        if (Check_Valid_Helix(residueList[position])): 
            # print ("Check Valid 3")
            # Converting "-" to "H"
            for j in range (start, end+1):
                List_Containing_Possible_Helices[j] = "H"

            left = dict_with_residuePosn_and_indices[position][0] - 1         # Index of the front pointer, to be moved backwards
            right = dict_with_residuePosn_and_indices[position][1] + 1           # Index of the last pointer, to be moved frontwards
        
        # Case 1 : Moving rightwards
            current_residue = residueList[position][-3:] + sequence_given[right]
            # Assuming residue to be helical, hence assigning.
            Helical_Region = residueList[position]
            while (right<len(sequence_given)):
                if (Helix_Score_using_P(current_residue) >= 4):
                    # print ("Check 8")
                    Helical_Region += current_residue[-1]
                    if (right+1 != len(sequence_given)):
                        List_Containing_Possible_Helices[right+1] = "H"
                        current_residue = current_residue[-3:] + sequence_given[right+1]
                    right+=1
                else :
                    break
        # Case 2 : Moving leftwards; updating currrent_residue, while Helical region updated waala hi use kareinge, naaki usse reset kar denge to residueList[position]
            current_residue =  sequence_given[left] + residueList[position][ : 3]
            while (left>=0):
                if (Helix_Score_using_P(current_residue) >= 4):
                    # print ("Check 9")
                    Helical_Region = current_residue[0] + Helical_Region
                    List_Containing_Possible_Helices[left] = "H"
                    left-=1
                    current_residue =  sequence_given[left] + current_residue[ : 3]
                else : 
                    break

            List_Of_Helical_Regions.append(Helical_Region)

# print(residueList)
# print(List_Of_Helical_Regions)
print("Considering cases only of Helical Regions : \n", "".join(List_Containing_Possible_Helices), end = "\n\n")



# Doing analysis for beta
def return_possibilities_beta (sequence) : 
    dict_with_residue_and_indices = {}
    possibilities = []
    for i in range(len(sequence)-4):
        possibilities.append(sequence[i : 5+i])
        dict_with_residue_and_indices[i] = [i, 5+i -1]         # i represents the nth index in the possibilities list.
    return possibilities, dict_with_residue_and_indices

# Returns score of the string entered for beta
def Beta_Score_using_P (residue) :    
    score = 0
    for aminoAcid in residue:
        for key in amino_acid_values :
            if (key == aminoAcid):
                score+= amino_acid_values[aminoAcid][1]
                # 1st index contains P(s)
    return score

def Check_Valid_Beta (residue):
    count = 0
    for aminoAcid in residue:
        for key in amino_acid_values :
            if (key == aminoAcid):
                if(amino_acid_values[aminoAcid][1] >= 1) :
                    count +=1
    if (count >= 3) : 
        return True
    return False

residueList_beta, dict_with_residuePosn_and_indices_beta = return_possibilities_beta(sequence_given)
# print(return_possibilities_beta(sequence_given))

List_Of_Beta_Regions = []
List_Containing_Possible_Betas = ["-" for i in range(len(sequence_given))]


for position in dict_with_residuePosn_and_indices_beta:
    start, end = position, dict_with_residuePosn_and_indices_beta[position][1] # Both included.
    # Case when first residue => neglecting the effect of front pointer.
    if (position == 0): 
        # print ("Check 1")
        # Primary check
        if (Check_Valid_Beta(residueList_beta[position])): 
            # print ("Check Valid 1")
            # Converting "-" to "H"
            for j in range (start, end+1):
                List_Containing_Possible_Betas[j] = "S"

            i=4+1
            current_residue = residueList_beta[position][-3:] + sequence_given[i]
            # Assuming residue to be beta one, hence assigning.
            Beta_Region = residueList_beta[position]
            while (i<len(sequence_given)):
                if (Beta_Score_using_P(current_residue) >= 4):
                    # print ("Check 2")
                    Beta_Region += current_residue[-1]
                    if (i+1 != len(sequence_given)):
                        List_Containing_Possible_Betas[i+1] = "S"
                        current_residue = current_residue[-3:] + sequence_given[i+1]
                    i+=1
                else :     
                    break
            List_Of_Beta_Regions.append(Beta_Region)

# Case last residue => neglecting the effect of end pointer.
    elif (position == len(residueList_beta)-1) :
        # print ("Check 4")
        # Primary check
        if (Check_Valid_Beta(residueList_beta[position])): 
            # print ("Check Valid 2")
            # Converting "-" to "H"
            for j in range (start, end+1):
                List_Containing_Possible_Betas[j] = "S"

            i = dict_with_residuePosn_and_indices_beta[len(residueList_beta)-1][0] - 1
            current_residue =  sequence_given[i] + residueList_beta[position][ : 3]
            # Assuming residue to be beta, hence assigning.
            Beta_Region = residueList_beta[position]
            while (i>=0):
                if (Beta_Score_using_P(current_residue) >= 4):
                    # print ("Check 5")
                    Beta_Region = current_residue[0] + Beta_Region
                    List_Containing_Possible_Betas[i] = "S"
                    i-=1
                    current_residue =  sequence_given[i] + current_residue[ : 3]
                else : 
                    break
            
            List_Of_Beta_Regions.append(Beta_Region)

    # Case of residues other than 1st and last.
    else :
        # print ("Check 7")
        # Primary check
        if (Check_Valid_Beta(residueList_beta[position])): 
            # print ("Check Valid 3")
            # Converting "-" to "H"
            for j in range (start, end+1):
                List_Containing_Possible_Betas[j] = "S"

            left = dict_with_residuePosn_and_indices_beta[position][0] - 1         # Index of the front pointer, to be moved backwards
            right = dict_with_residuePosn_and_indices_beta[position][1] + 1           # Index of the last pointer, to be moved frontwards
        
        # Case 1 : Moving rightwards
            current_residue = residueList_beta[position][-3:] + sequence_given[right]
            # Assuming residue to be beta, hence assigning.
            Beta_Region = residueList_beta[position]

            while (right<len(sequence_given)):
                if (Beta_Score_using_P(current_residue) >= 4):
                    # print ("Check 8")
                    Beta_Region += current_residue[-1]
                    if (right+1 != len(sequence_given)):
                        List_Containing_Possible_Betas[right+1] = "S"
                        current_residue = current_residue[-3:] + sequence_given[right+1]
                    right+=1
                else :
                    break
        # Case 2 : Moving leftwards; updating currrent_residue, while beta region updated waala hi use kareinge, naaki usse reset kar denge to residueList_beta[position]
            current_residue =  sequence_given[left] + residueList_beta[position][ : 3]
            while (left>=0):
                if (Beta_Score_using_P(current_residue) >= 4):
                    # print ("Check 9")
                    Beta_Region = current_residue[0] + Beta_Region
                    List_Containing_Possible_Betas[left] = "S"
                    left-=1
                    current_residue =  sequence_given[left] + current_residue[ : 3]
                else : 
                    break

            List_Of_Beta_Regions.append(Beta_Region)

# print(residueList_beta)
# print(List_Of_Beta_Regions)
print("Considering cases only of Beta Regions : \n", "".join(List_Containing_Possible_Betas), end = "\n\n")
# print("List of beta regions :", List_Of_Beta_Regions)



# ------------------------------------------------------------------

# Resolving Conflicts :
# Making a dictionary, which contains key as the index, and value as the list of charachters on that index in the 2 lists

dict_conflicts_resolve = {}
Conflicting_Regions = []
resolved_List = [0]*len(sequence_given)

for i in range(len(sequence_given)):
    dict_conflicts_resolve[i] = [List_Containing_Possible_Helices[i], List_Containing_Possible_Betas[i]]

i = 0
while (i<len(sequence_given)): 
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
        
print("Portions of Conflict : \n", Conflicting_Regions, end = "\n\n")
print("Final Answer after Resolving Conflicts : \n", "".join(resolved_List), end = "\n\n")
