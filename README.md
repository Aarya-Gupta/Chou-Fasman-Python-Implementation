# Chou Fasman Protein Secondary Structure Prediction - Python Implementation

This code implements the Chou Fasman method to predict the secondary structure (helices and strands) of a protein sequence.

**Explanation (only for code.py):**

1. **Data and Initialization:**
    * A dictionary `amino_acid_values` stores scores for each amino acid related to helix (P(h)) and strand (P(s)) propensity.
    * A protein sequence `sequence_given` is defined.
    * Three lists are initialized:
        * `List_Containing_Possible_Helices`: Stores predicted helix residues ('-'' for non-helix).
        * `List_Containing_Possible_Betas`: Stores predicted strand residues ('-'' for non-strand).
        * `output`: Stores the final predicted secondary structure ('-'' for undetermined).

2. **Scoring Functions:**
    * `Helix_Score_using_P(residue)`: Calculates the total helix propensity score for a given amino acid sequence `residue`.
    * `Beta_Score_using_P(residue)`: Calculates the total strand propensity score for a given amino acid sequence `residue`.

3. **Nucleation Site Calculation:**
    * `calculate(i, type)`: Analyzes a window of amino acids around position `i` to determine if it's a valid nucleation site for a helix or strand based on the `type` ("S" for strand, "H" for helix).

4. **Extension Functions:**
    * `extend_left(j, type)`: Extends a predicted helix/strand region (`type`) to the left from position `j` as long as the P(h)/P(s) score remains favorable.
    * `extend_right(j, type)`: Extends a predicted helix/strand region (`type`) to the right from position `j` as long as the P(h)/P(s) score remains favorable.

5. **Sequence Processing:**
    * `process_sequence(sequence, target, length)`: Processes the entire sequence recursively. It calls `calculate` to identify potential nucleation sites for the target secondary structure (`target`, "H" for helix, "S" for strand) with a minimum length of `length`. If a valid site is found, it extends the predicted region using `extend_left` and `extend_right`.

6. **Conflict Resolution:**
    * A dictionary `dict_conflicts_resolve` is created to store predicted helix and strand characters at each position.
    * A list `Conflicting_Regions` stores the start and end indices of regions with conflicting helix and strand predictions.
    * A list `resolved_List` is created to store the final predicted secondary structure for each amino acid.
    * The code iterates through the sequence and resolves conflicts based on the following logic:
        * If neither helix nor strand is predicted, the residue is marked as '-'.
        * If only one (helix or strand) is predicted, the residue is marked accordingly.
        * If both helix and strand are predicted for a region, the amino acid sequence within that region is evaluated using the scoring functions `Helix_Score_using_P` and `Beta_Score_using_P`. The structure with the higher score is assigned to the conflicting region in `resolved_List`. If the scores are equal, the code prioritizes strands by default.

7. **Output:**
    * The code prints the following information:
        * Predicted helices only.
        * Predicted strands only.
        * List of conflicting regions.
        * Final predicted secondary structure after resolving conflicts.

**Overall, this code implements the Chou Fasman method to predict the secondary structure of a protein sequence. It considers the amino acid propensities for helix and strand formation and employs a window-based approach to identify potential secondary structure elements.**

**Note:** The Chou Fasman method is a relatively simple approach with limitations. More sophisticated methods might offer higher accuracy for protein secondary structure prediction.
