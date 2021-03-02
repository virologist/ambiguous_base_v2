# _*_ coding: utf-8 _*_
import pandas as pd
import numpy as np
import os
import sys
from Bio import SeqIO
from collections import Counter
import os

print(
'''
Date: 2021-1-28
Author: Yusy
Version: V0.2
Usage: python ambiguous_base.py filename
''')
# 删除clean_data.fasta文件
if os.path.exists("clean_data.fasta"):
    os.remove("clean_data.fasta")

# 计算汉明距离
def hammingDistance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

# 读取序列
print("=====================正在读取序列=====================")
filename = sys.argv[1]
Sequences = {}
Names = []
for seq_record in SeqIO.parse(filename, "fasta"):
    Sequences[seq_record.id] = str(seq_record.seq).upper()
    Names.append(seq_record.id)
print("=====================序列读取完成=====================")

Seqlen = len(Sequences[Names[0]])
for i in Sequences.values():
    assert len(i) == Seqlen , "序列没有经过比对裁剪"

print("=====================正在创建汉明距离矩阵（这一步需要较长时间）=====================")
# 构建汉明距离矩阵
df = pd.DataFrame(columns=Names,index=Names)
for i in range(len(df.columns)):
    for j in range(len(df.index)):
        df.iloc[i,j] = hammingDistance(Sequences[df.columns[i]],Sequences[df.index[j]])
df.to_csv('hamming_distance.csv')

print("=====================创建汉明距离矩阵完成=====================")
# 读取汉明距离矩阵
df = pd.read_csv('hamming_distance.csv', index_col=0)
# print(df.head())

# 模糊碱基及含义
ambiguous_bases = {'R': ['A', 'G'],'Y': ['C', 'T', 'U'],'M': ['A', 'C'],'K': ['G', 'T'],'S': ['C', 'G'],'W': ['A', 'T'],'H': ['A', 'C', 'T'],'B': ['C', 'G', 'T'],'V': ['A', 'C', 'G'],'D': ['A', 'G', 'T'],'N': ['A', 'C', 'G', 'T']}

# 逐条序列寻找模糊碱基
get_ab = {}
for name in Names:
    i = 0
    ab_position_list = []
    ab_name_list = []
    for i in range(len(Sequences[name])):
        if Sequences[name][i] in ambiguous_bases.keys():
            print("{}这条序列的第{}号位点有模糊碱基{}".format(name,str(i+1),Sequences[name][i]))
            ab_position_list.append(str(i))
            ab_name_list.append(Sequences[name][i])
    get_ab[name] = ab_position_list
    i += 1
# 验证匹配到的数目和位点是否正确
# test = Counter(Sequences['KF226105'])
# print(test)
# print(Sequences['KF226105'].index('W'))
get_ab_not_empty = {}
for key,value in get_ab.items():
    if len(value) != 0:
        get_ab_not_empty[key] = value
#print("模糊碱基包括:{}".format(get_ab_not_empty))
for name,position in get_ab_not_empty.items():
    count = 1
    print("{:=^50s}".format("Split Line"))
    print("{}共有{}个模糊碱基，位点是:{}".format(name,len(position),position))
    top_10 = df.sort_values(by=name).head(11)[name]
    top_10 = top_10[1:10]
    top_10 = list(top_10.index)#获取到汉明距离最近的十个序列的ID
    #print("与{}汉明距离最近的十条序列分别是：{}".format(name,top_10))
    for i in position:
        #print(i)
        position_base = []
        for j in top_10:
            position_base.append(Sequences[j][int(i)])
        #("汉明距离最近的十条序列该位点的碱基分别是：{}".format(position_base))
        position_base_count = Counter(position_base)
        if position_base_count.most_common(1)[0][0] in ambiguous_bases[Sequences[name][int(i)]]:
            base_replace = position_base_count.most_common(1)[0][0]
        else:
            base_replace = ambiguous_bases[Sequences[name][int(i)]][0]
        print("*第{}位的碱基{}将被替换成碱基{}".format(i,Sequences[name][int(i)],base_replace))
        #为了解决str不可更改，引入中间变量（列表）
        Med_list = []
        Med_list = list(Sequences[name])
        Med_list[int(i)] = base_replace
        Sequences[name] = "".join(Med_list)

#清洗提前终止密码子
st_codon_to = []
all_words = []
st_codon_position = []
for value in Sequences.values():
    plus = str(value)
    words = [plus[s:s+3] for s in range(0,len(value)-3,3)]
    all_words.append(words)

o = 0
for i in range(len(all_words)):
    for j in range(len(all_words[i])):
        if all_words[i][j] in ["TAG","TGA","TAA"]:
            o += 1
            st_codon_position.append(str(i))
            print("=====================发现提前终止密码子，正在清洗=====================")
            if i < 3:
                st_codon_to.append(all_words[i + 1][j])
                st_codon_to.append(all_words[i + 2][j])
                st_codon_to.append(all_words[i + 3][j])
                st_codon_to.append(all_words[i + 4][j])
                st_codon_to_count = Counter(st_codon_to)
                top_codon = st_codon_to_count.most_common(1)[0][0]
                all_words[i][j] = top_codon
            else:
                st_codon_to.append(all_words[i - 1][j])
                st_codon_to.append(all_words[i - 2][j])
                st_codon_to.append(all_words[i + 1][j])
                st_codon_to.append(all_words[i + 2][j])
                st_codon_to_count = Counter(st_codon_to)
                top_codon = st_codon_to_count.most_common(1)[0][0]
                all_words[i][j] = top_codon
            name = Names[i]
            Sequences[name] = "".join(all_words[i])
print("共有{}个提前终止的密码子，出现在第{}条序列".format(str(o),st_codon_position))

with open('clean_data.fasta', 'a') as f:
    for key,value in Sequences.items():
        f.write('>'+key)
        f.write('\n')
        f.write(value)
        f.write('\n')
print("=====================Done=====================")

