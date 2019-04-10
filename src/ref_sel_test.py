import codecs
import csv
import pandas
import sys, operator


def getRefExpLevel(lowerTh, upperTh, outputfile):
    with codecs.open("final_out.csv", "r", encoding='utf-8', errors='ignore') as f:
        with open(outputfile, 'w') as csvfile:
          reader = csv.DictReader(f)
          writer = csv.writer(csvfile)
          writer.writerow(["genes", "cov", "mean", "std", "MFC", "covMFC"])
          for row in reader:
            #try:

                ens = []
                if lowerTh <= pandas.to_numeric(row["mean"]) <= upperTh:
                  covMFC = float(row["MFC"])*float(row["cov"])
                  Account = [row["genes"], row["cov"], row["mean"], row["std"], row["MFC"], covMFC]
                #print(row["Gene_ID"] + " " + mg.getgene(row["Gene_ID"])['symbol'])
                  writer.writerow(Account)

def findTop(resultFile):
    df  = pandas.read_csv(resultFile)
    sortedlist = df.sort_values(by='covMFC')
    for index, row in sortedlist.iterrows():
        if row[0].rsplit('.',1)[0] in lines:
            print(row[0])
            break


text_file = open("data/house_keeping_ensembl.txt", "r")
lines = text_file.read().splitlines()

outfiles = ['ref1.csv','ref2.csv','ref3.csv','ref4.csv','ref5.csv']
lows = [3,14,50,220,800]
ups = [5,18,80,280,1200]

for i in range(len(outfiles)):
    getRefExpLevel(lows[i],ups[i],outfiles[i])

for outfile in outfiles:
    findTop(outfile)

sys.exit(0)
