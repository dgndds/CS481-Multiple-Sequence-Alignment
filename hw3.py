from nw import needleman_wunsch



def scoreColumnAlignment(char,profile,j,match,mismatch,gap):
    score = 0
    alphabet = "ACTG-"

    for i in range(0,5):
        # if char == alphabet[i] and char != "-" and alphabet[i] != "-":
        #     score = score + (match * getCharFreq(alphabet[i],profile,j))
        # elif char !=alphabet[i]:
        #     if char != "-"  and alphabet[i] != "-":
        #         score = score + (mismatch * getCharFreq(alphabet[i],profile,j))
        #     elif char == "-" or alphabet[i] == "-":
        #         score = score + (gap * getCharFreq(alphabet[i],profile,j))
        if char == alphabet[i] and char != "-":
            score = score + (match * getCharFreq(alphabet[i],profile,j))
        elif char !=alphabet[i] or char == "-":
            if char != "-":
                score = score + (mismatch * getCharFreq(alphabet[i],profile,j))
            elif char == "-":
                score = score + (gap * getCharFreq(alphabet[i],profile,j))    
        
    return score

def generateProfile(seqs):
    m = len(seqs[0])
    print(m)
    scoreTable = [[0 for i in range(m)] for j in range(5)] 

    for i in range(0,m):
        counts = [0 for i in range(5)]
        
        for j in range (0,5):
            if(seqs[j][i] == "A"):
                counts[0] = counts[0] + 1
            elif(seqs[j][i] == "C"):
                counts[1] = counts[1] + 1
            elif(seqs[j][i] == "T"):
                counts[2] = counts[2] + 1
            elif(seqs[j][i] == "G"):
                counts[3] = counts[3] + 1
            elif(seqs[j][i] == "-"):
                counts[4] = counts[4] + 1

        for j in range(0,5):
            scoreTable[j][i] = counts[j] / 5
    

    return scoreTable

def getCharFreq(char,profile,index):
    charIndex = "ACTG-".index(char)
    return profile[charIndex][index]

def sequenceToProfileAlign(seq1, seq2, match=2, mismatch=-1, gap=-1):
    # seq1 = seqs
    # seq2 = seq
    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -1

    m, n = len(seq1[0]), len(seq2)

    V = [[0 for i in range(n+1)] for j in range(m+1)]
    P = [[0 for i in range(n+1)] for j in range(m+1)]

    for i in range(m+1): # seqs

        for j in range(n+1): #seq

            if i >= 0 and j == 0:
                V[i][j] = gap_penalty * i
                P[i][j] = 1

            elif i == 0 and j >= 0:
                V[i][j] = gap_penalty * j
                P[i][j] = 2

            else:

                # if seq1[i-1] == seq2[j-1]:
                match_case = V[i-1][j-1] + scoreColumnAlignment(seq2[j-1],profile,i-1,match,mismatch,gap)
                # else:
                    # match_case = V[i-1][j-1] + mismatch_penalty

                gap_in_first_case = V[i][j-1] + gap_penalty
                

                gap_in_second_case = V[i-1][j] + scoreColumnAlignment("-",profile,i-1,match,mismatch,gap)

                value_list = [match_case,
                              gap_in_first_case, gap_in_second_case]
                max_value = value_list.index(max(value_list))
                V[i][j] = value_list[max_value]
                P[i][j] = max_value

    print(V)
    # backtracking from the last indices
    # i = m
    # j = n

    # aligned_seq1 = ""
    # aligned_seq2 = ""

    # score = 0

    # while i > 0 and j > 0:

    #     # match or mismatch
    #     if P[i][j] == 0:
    #         aligned_seq1 += seq1[i-1]
    #         aligned_seq2 += seq2[j-1]
    #         i -= 1
    #         j -= 1
    #         score += match_score

    #     # gap in the second seq
    #     elif P[i][j] == 2:
    #         aligned_seq1 += seq1[i-1]
    #         aligned_seq2 += "-"
    #         i -= 1
    #         score += gap_penalty

    #     # gap in the first seq
    #     elif P[i][j] == 1:
    #         aligned_seq1 += "-"
    #         aligned_seq2 += seq2[j-1]
    #         j -= 1
    #         score += gap_penalty

    # if i > 0:
    #     aligned_seq1 += seq1[i:]
    #     aligned_seq2 += "-"
    #     i -= 1
    #     score += gap_penalty

    # if j > 0:
    #     aligned_seq1 += "-"
    #     aligned_seq2 += seq2[j:]
    #     j -= 1
    #     score += gap_penalty

    # # read from the end
    # seq1 = ''.join(reversed(aligned_seq1))
    # seq2 = ''.join(reversed(aligned_seq2))

    # return [score, seq1, seq2]

seqs = [
    "-AGGCTATCACCTG",
    "TAG-CTACCA---G",
    "CAG-CTACCA---G", 
	"CAG-CTATCAC-GG", 
	"CAG-CTATCGC-GG" 
]

seq = "CAGGTACCACGG"

# print(seqs[2][0])

profile = generateProfile(seqs)
# print()
# print(getCharFreq("C",profile,0))
# print(scoreColumnAlignment("C",profile,0,2,-1,-1))
sequenceToProfileAlign(seqs,seq,2,-1,-1)