# libiary
import numpy as np
import pandas as pd


def alignment(input_path, score_path, output_path, aln, gap):

    # scoreTable
    pam = pd.read_csv(score_path, sep=r'\s+', skiprows=9)
    oddmatrix = np.array(pam)
    
    # readFile
    name = []  # sequencesName
    data = []  # sequenceData
    inFIle = pd.read_csv(input_path, header=None)
    for index, row in inFIle.iterrows():
        if index % 2 == 0:
            name.append(row[0])
        else: 
            data.append(row[0])
    index += 1 # totalLine

    # global or local
    if aln == "global":
        globalPoint(name, data, oddmatrix, gap, output_path)
    else:
        localPoint(name, data, oddmatrix, gap, output_path)
    


    return 0

def globalPoint(names, data, score, gap, output_path):
    # create point and step matrix
    point = np.zeros((len(data[0])+1, len(data[1])+1))
    step = np.zeros((len(data[0])+1, len(data[1])+1))
    # init
    for i in range(1, len(data[0])+1):
        point[i][0] = point[i-1][0] + gap
        step[i][0] = 3
    for j in range(1, len(data[1])+1):
        point[0][j] = point[0][j-1] + gap
        step[0][j] = 2

    # use dp count score, step tatrix : left 2, right 3, match 5

    for i in range(1, len(data[0])+1):
        for j in range(1, len(data[1])+1):
            if point[i-1][j] > point[i][j-1]:
                if point[i-1][j] + gap > point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]:
                    point[i][j] = point[i-1][j] + gap
                    step[i][j] = 3
                else:
                    point[i][j] = point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]
                    step[i][j] = 5
            # elif (point[i-1][j] == point[i][j-1]) and point[i-1][j] + gap > point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]:
            #     point[i][j] = point[i][j-1] + gap
            #     step[i][j] = 2 * 3
            else:
                if point[i][j-1] + gap > point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]:
                    point[i][j] = point[i][j-1] + gap
                    step[i][j] = 2
                else:
                    point[i][j] = point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]
                    step[i][j] = 5
    
    # collect sequence
    seq1 = ""
    seq2 = ""
    while ((i > 0) or (j > 0)):
        if step[i][j] == 2:
            seq1 += '-'
            seq2 += data[1][j-1]
            j -= 1
        elif step[i][j] == 3:
            seq1 += data[0][i-1]
            seq2 += '-'
            i -= 1
        elif step[i][j] == 5:
            seq1 += data[0][i-1]
            seq2 += data[1][j-1]
            i -= 1
            j -= 1

    seq1 = seq1[::-1]
    seq2 = seq2[::-1]
    
    # output file
    lines = [names[0], seq1, names[1], seq2]

    with open(output_path, "w") as file:
        for line in lines:
            file.write(line + "\n")
    return 0

def localPoint(names, data, score, gap, output_path):
    # create point and step matrix
    point = np.zeros((len(data[0])+1, len(data[1])+1))
    step = np.zeros((len(data[0])+1, len(data[1])+1))
    max = 0
    count = 0
    index = []

    # use dp count score, step tatrix : left 2, right 3, match 5, left and right 6, stop 0

    for i in range(1, len(data[0])+1):
        for j in range(1, len(data[1])+1):
            if point[i-1][j] > point[i][j-1]:
                if point[i-1][j] + gap > point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]:
                    if point[i-1][j] + gap >= 0:
                        point[i][j] = point[i-1][j] + gap
                        step[i][j] = 3
                        if point[i][j] > max:
                            max = point[i][j]
                            index = []
                            index.append([i, j])
                            count = 1
                        elif point[i][j] == max:
                            index.append([i, j])
                            count += 1
                    else:
                        point[i][j] = 0
                        step[i][j] = 0
                else:
                    if point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]] >= 0:
                        point[i][j] = point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]
                        step[i][j] = 5
                        if point[i][j] > max:
                            max = point[i][j]
                            index = []
                            index.append([i, j])
                            count = 1
                        elif point[i][j] == max:
                            index.append([i, j])
                            count += 1
                    else:
                        point[i][j] = 0
                        step[i][j] = 0
            elif (point[i-1][j] == point[i][j-1]) and point[i-1][j] + gap > point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]:
                if point[i-1][j] + gap >= 0:
                    point[i][j] = point[i-1][j] + gap
                    step[i][j] = 2 * 3
                    if point[i][j] > max:
                        max = point[i][j]
                        index = []
                        index.append([i, j])
                        count = 1
                    elif point[i][j] == max:
                        index.append([i, j])
                        count += 1
                else:
                    point[i][j] = 0
                    step[i][j] = 0
            else:
                if point[i][j-1] + gap > point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]:
                    if point[i][j-1] + gap >= 0:
                        point[i][j] = point[i][j-1] + gap
                        step[i][j] = 2
                        if point[i][j] > max:
                            max = point[i][j]
                            index = []
                            index.append([i, j])
                            count = 1
                        elif point[i][j] == max:
                            index.append([i, j])
                            count += 1
                    else:
                        point[i][j] = 0
                        step[i][j] = 0
                else:
                    point[i][j] = point[i-1][j-1] + score[np.where(aminoAcid == data[0][i-1])[0][0]][np.where(aminoAcid == data[1][j-1])[0][0]]
                    step[i][j] = 5
                    if point[i][j] > max:
                        max = point[i][j]
                        index = []
                        index.append([i, j])
                        count = 1
                    elif point[i][j] == max:
                        index.append([i, j])
                        count += 1
                    if point[i][j] < 0:
                        point[i][j] = 0
                        step[i][j] = 0
    
    # find max score and longest sequences
    ansName = []
    print(index)
    for k in index:
        seq1 = ""
        seq2 = ""
        ansName = findLongPair(ansName, seq1, seq2, k[0], k[1], step, data)
    max = 0
    for k in ansName:
        if len(k[0]) > max:
            max = len(k[0])

    # output file
    for k in ansName:
        if len(k[0]) == max:
            k[0] = k[0][::-1]
            k[1] = k[1][::-1]
            lines = [names[0], k[0], names[1], k[1]]
            with open(output_path, "w") as file:
                for line in lines:
                    file.write(line + "\n")
    
    
    return 0

# collect all max score and longest sequences
def findLongPair(ansName, seq1, seq2, i, j, step, data):
    while ((i > 0) or (j > 0)):
        if step[i][j] == 2:
            seq1 += '-'
            seq2 += data[1][j-1]
            j -= 1
        elif step[i][j] == 3:
            seq1 += data[0][i-1]
            seq2 += '-'
            i -= 1
        elif step[i][j] == 5:
            seq1 += data[0][i-1]
            seq2 += data[1][j-1]
            i -= 1
            j -= 1
        elif step[i][j] == 6:
            seq1 += '-'
            seq2 += data[1][j-1]
            ansName = findLongPair(seq1, seq2, i, j-1, step, data)
            seq1 = seq1[:-1]
            seq2 = seq2[:-1]
            seq1 += data[0][i-1]
            seq2 += '-'
            ansName = findLongPair(seq1, seq2, i-1, j, step, data)
        else:
            break
    
    ansName.append([seq1, seq2])
    return ansName

#main

aminoAcid = np.array(["A", "R", "N", "D", "C",
                      "Q", "E", "G", "H", "I",
                      "L", "K", "M", "F", "P",
                      "S", "T", "W", "Y", "V"])

# alignment("BioInfo/HW3/examples/test_global.fasta","BioInfo/HW3/examples/pam250.txt","BioInfo/HW3/result.fasta","global",-10)
# alignment("BioInfo/HW3/examples/test_local.fasta","BioInfo/HW3/examples/pam100.txt","BioInfo/HW3/result.fasta","local",-10)