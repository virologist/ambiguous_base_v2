##################使用须知#####################
1.单条序列中歧义碱基的含量不应大于1%或2%，不推荐歧义碱基含量高的序列使用
2.校正歧义碱基的算法是根据汉明距离最近的序列，因此建议下载所有可得序列后进行去歧义碱基计算，最后再考虑去除重复或者相似度高的序列。
3.此脚本不合适校正遗传距离较远的序列歧义碱基位点
4.输入序列需经过比对和裁剪保留ORF，第一位密码子需位ATG三联体密码子。
5.在满足4的条件下，异常终止密码子为提前终止的三联体密码子，通常是由测序质量不佳，在拼接过程中出现的。

######
Usage: python ambiguous_base.py input_align_sequences.fas
Becareful, the input sequences (fas/fasta format) shoud be align and trim
Function: correct the ambiguous base

