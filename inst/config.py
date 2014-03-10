#!/usr/bin/python
#import global_file
import sys
import getopt
import csv
import argparse
import os, errno

def readValues(argv):
    #First getting the nPops
    parser = argparse.ArgumentParser(description='This is a demo script')
    parser.add_argument('-n','--ne_int', help='Ne', required=False)
    parser.add_argument('-o','--ne_mat', help='Ne', required=False)
    parser.add_argument('-p', '--pops', help='Starting Population Matrix', required=True)
    parser.add_argument('--rMean_mat', help='rMean', required=False)
    parser.add_argument('--rMean_int', help='rMean', required=False)
    parser.add_argument('--rVar_mat', help='rVar', required=False)
    parser.add_argument('--rVar_int', help='rVar', required=False)
    parser.add_argument('--K_int', help='Kint', required=False)
    parser.add_argument('--K_mat', help='Kmat', required=False)
    parser.add_argument('--A_int', help='Aint', required=False)
    parser.add_argument('--A_mat', help='Amat', required=False)
    parser.add_argument('--catProb_int', help='catProb', required=False)
    parser.add_argument('--catProb_mat', help='catProb', required=False)
    parser.add_argument('--diploid', help='diploid or haploid', required=True)
    parser.add_argument('--nGens', help='nGens', required=True)
    parser.add_argument('--mig_mat', help='mig mat', required=False)
    parser.add_argument('--mig_int', help='mig int', required=False)
    parser.add_argument('--SNP_model_input', help='SNP_model_input', required=True)
    parser.add_argument('--h_input', help='h_input', required=True)
    parser.add_argument('--s_input', help='s_input', required=False)
    parser.add_argument('--s_mat', help='s_mat', required=False)
    parser.add_argument('--gamma_alpha_input', help='gamma_alpha_input', required=False)
    parser.add_argument('--gamma_beta_input', help='gamma_beta_input', required=False)
    parser.add_argument('--gSize_input', help='gSize', required=False)
    parser.add_argument('--mutRate_input', help='mutRate', required=False)
    parser.add_argument('--GENEPOP', help='GENEPOP', required=False)
    parser.add_argument('--GENELAND', help='GENELAND', required=False)
    parser.add_argument('--PLINK', help='PLINK', required=False)
    parser.add_argument('--outfile', help='outfile', required=True)
    parser.add_argument('--nSNPs_input', help='outfile', required=False)
    parser.add_argument('--SNPs_starting_freq', help='outfile', required=False)
    parser.add_argument('--infile', help='infile from previous run', required=False)
    parser.add_argument('--recordTrag', help='Whether or not you want to record allele tragectories', required=True)
    parser.add_argument('--sDiff', help='Matrix of s values, if they are different in different populations', required=False)
    
    args = parser.parse_args()
    
    mat = readFile(args.pops, remove=True)
    makePopulations(mat)
    
    global recordTrag
    recordTrag = int(args.recordTrag)

    global sDiff
    sDiff = None
    if args.sDiff != None:
        sDiff = readFile(args.sDiff, remove=True)
        #IMPLEMENT SOMETHING TO CHECK THAT THEY DON'T INPUT S and SDIFF

    ##Fix the incompatiblities
    global inFile_YEP
    inFile_YEP = None
    if args.infile != None:
        readInputfile(args.infile)
    
    global SNP_model
    SNP_model = int(args.SNP_model_input) #1='snm', 0='fixed'
    if SNP_model != 1 and SNP_model !=0:
        print("SNP MODEL MUST BE 0 or 1")
        sys.exit()
    
    global SNPs_starting_freq
    if args.SNPs_starting_freq != None and args.infile == None:
        SNPs_starting_freq= float(args.SNPs_starting_freq)
    
    global nSNPs
    if args.nSNPs_input != None and args.infile == None:
        nSNPs = int(args.nSNPs_input)
    
    global mutRate
    if args.mutRate_input != None:
        mutRate = float(args.mutRate_input)
    else:
        mutRate = None

    global gSize
    if args.gSize_input != None:
        gSize = int(float(args.gSize_input))
    else:
        gSize = None
    
    global GENELAND
    if args.GENELAND == None or args.GENELAND == "F" or args.GENELAND == "FALSE":
        GENELAND = False
    else:
        GENELAND = True

    global GENEPOP
    if args.GENEPOP == None or args.GENEPOP == "F" or args.GENEPOP == "FALSE":
        GENEPOP = False
    else:
        GENEPOP = True

    global PLINK
    if args.PLINK == None or args.PLINK == "FALSE" or args.PLINK == "F":
        PLINK = False
    else:
        PLINK = True
    
    global outFile
    outFile = args.outfile
    
    
    global nGens
    nGens = int(args.nGens)
    
    global h
    h = float(args.h_input)

    global s_mat, s, g_a, g_b
    if args.s_input != None:
        s_mat = None
        s = float(args.s_input)
    elif args.s_mat != None:
        s_mat = readFile(args.s_mat, remove=True)
        s = None
    else:
        g_a = float(args.gamma_alpha_input)
        g_b = float(args.gamma_beta_input)
        gSize = int(args.gSize)
        mutRate = float(args.mutRate)

    global migProb, mig #MIGRATION
    if args.mig_mat != None:
        migProb = readFile(args.mig_mat, remove=True)
        migProb = [map(float,x) for x in migProb]
        mig = True
    else:
        if args.mig_int == 0:
            mig = False
            migProb = 0
        else:
            mig = True
            migProb = float(args.mig_int)

    global Ne
    if args.ne_mat != None:
        Ne = readFile(args.ne_mat, remove=True)
        Ne = [map(int,x) for x in Ne]
    elif args.ne_int != None:
        Ne = int(args.ne_int)
    #else:
    #print("Must have an Ne input")
    #sys.exit()
    global diploid
    if args.diploid == False:
        diploid = False
    else:
        diploid = True
    

    global catProb
    if args.catProb_mat != None:
        catProb = readFile(args.catProb_mat, remove=True)
        catProb = [map(float,x) for x in catProb]
    elif args.catProb_int != None:
        catProb = int(args.catProb_int)

    global rMean, rVar
    if args.rMean_mat != None:
        rMean = readFile(args.rMean_mat, remove=True)
        rMean = [map(float,x) for x in rMean]
    elif args.rMean_int != None:
        rMean = float(args.rMean_int)
    
    if args.rVar_mat != None:
        rVar = readFile(args.rVar_mat, remove=True)
        rVar = [map(float,x) for x in rVar]
    
    elif args.rVar_int != None:
        rVar = float(args.rVar_int)

    global A, K
    if args.A_mat != None:
        A = readFile(args.A_mat, remove=True)
        A = [map(int,x) for x in A]
    elif args.A_int != None:
        A = int(args.A_int)
    
    if args.K_mat != None:
        K = readFile(args.K_mat, remove=True)
        K = [map(int,x) for x in K]
    elif args.K_int != None:
        K = int(args.K_int)
    
    global allele_history
    allele_history = list()




def readFile( fileToRead, remove ):
    f = open( fileToRead, 'r')
    try:
        l = [ map(str , line.rstrip('\n').split(',')) for line in f ]
    except:
        print('Problem with pops CSV file!')
        sys.exit()
    finally:
        f.close()
    if remove == True:
        silentremove( fileToRead )
    return l

def makePopulations(mat):
    global max_X, max_Y, populations, popMat
    ##Now check format + make 2 sets of tuples for initial + final pops
    populations = list()
    #max_X = len(mat[0])
    #max_Y = len(mat)
    max_X = len(mat) #This is the number of rows.
    max_Y = len(mat[0]) #This is the number of columns
    popMat = mat
    if max_X == 1 and max_Y == 1:
        if int(mat[0][0]) == 1: #It is a currently existing population!
            populations.append( (0,0) )
        elif int(mat[0][0]) == 0: #It is not currently existing, but it could!
            populations.append( (0,0) )
        elif int(mat[0][0]) == -1: #It is a barrier and cannot be entereddd #This doesn't make sense, take away this option
            pass
        else:
            print("Error in pop matrix values.  Must be (1,0,-1)")
            sys.exit()
    else:
        # print str(len(mat[0])) + " " + str(len(mat))
        for j in range(0,len(mat[0])): #for each COLUMN
            for i in range(0,len(mat)): #for each ROW
                #print str(i) + " " + str(j)
                if int(mat[i][j]) == 1: #It is a currently existing population!
                    populations.append( (i,j) )
                elif int(mat[i][j]) == 0: #It is not currently existing, but it could!
                    populations.append( (i,j) )
                elif int(mat[i][j]) == -1: #It is a barrier and cannot be entereddd
                    populations.append( (i,j) ) #basically I'm not using the -1 for anything
                else:
                    print("Error in pop matrix values.  Must be (1,0,-1)")
                    sys.exit()

def readInputfile(infile):
    #print("AHOY")
    global nSNPs, Ne, popSizes, populations, max_X, max_Y, inFile_YEP
    # inFile_YEP = True
    inFile_YEP = infile
    
    f = open( infile, 'r')
    try:
        l = [ map(str , line.rstrip('\n').split(',')) for line in f ]
    except:
        print('Problem with pops CSV file!')
        sys.exit()
    finally:
        f.close()
    
    parsedFile = l
    #    max_X, max_Y, nSNPs = parsedFile[0][0:3] #This should be the max_X and max_Y
    
    max_X = int(parsedFile[0][0])
    max_Y = int(parsedFile[0][1])
    nSNPs = int(parsedFile[0][2])


##closed vs. open boundaries
def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured

