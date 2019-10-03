#!/usr/bin/python
########################################################################
"""
    <G4Hunter - a program to search quadruplex-forming regions in DNA.>
    Copyright (C) <2012>  <Bedrat amina  supervised by Dr.Jean Louis Mergny>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
########################################################################
import os, re, sys, getopt
import time
import shutil
import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
from Bio import SeqIO

 
def main(argv):

   if not argv:
        sys.stdout.write("Sorry: you must specify at least an argument\n")
        sys.stdout.write("More help avalaible with -h or --help option\n")
        sys.exit(1)

   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:w:s:",["help","ifile=","ofile="])
   except getopt.GetoptError:
      print '\033[1m' +'python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>\n'+'\033[0;0m'
      sys.exit(1)
      
   for opt, arg in opts:
      if opt in ('-h',"--help"):
          print  '\033[1m' +'\t ----------------------'+'\033[0;0m'
          print  '\033[1m' +'\n\t  Welcome To G4Hunter :'+'\033[0;0m'
          print  '\033[1m' +'\t ----------------------'+'\033[0;0m'

          print '\n G4Hunter takes into account G-richness and G-skewness of a given sequence and gives a quadruplex propensity score as output.'
          print 'To run G4Hunter use the commande line: \n'
          print  '\033[1m' +'python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>\n'+'\033[0;0m'
          sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile= arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-w", "--wind"):
         window = arg
      elif opt in ("-s", "--score"):
         score = arg
             
   return  inputfile, outputfile, int(window), float(score)

#calcule le score de chaque base dans un fichier qui peut contenir plusieur sequences fasta
#le fichier doit comporter la nom de la seq + la sequence en une seul ligne donc pas de \n ou \r 
class Soft(object):

    def __init__(self):
        pass
     
    def ReadFile(self,Filein):
        ListSeq,LHeader=[],[]
        for record in SeqIO.parse(Filein, "fasta") :
            LHeader.append(record.id)
            ListSeq.append(record.seq)
        return LHeader,ListSeq

    def ReadSeq(self, Seqs):
        ListSeq,LHeader=[],[]
        for record in SeqIO.parse(Seqs, "fasta") :
            LHeader.append(record.id)
            ListSeq.append(record.seq)
        return LHeader,ListSeq

    def GFinder(self,Filein,k):
        LHeader,ListSeq=self.ReadFile(Filein)
        LSeq,LNumber,LScoreSeq,SeqLine=[],[],[],""
        for i in range(len(ListSeq)) :
            Sequence,liste=self.BaseScore(ListSeq[i])
            LSeq.append(Sequence)
            LScoreSeq.append(self.CalScore(liste, k))
            LNumber.append(liste)
        return LScoreSeq, LSeq , LNumber, LHeader
    
    def BaseScore(self,line):
        item, liste=0, []
        #calcule le item de chaque base et la stock dans liste
        while ( item < len(line)):
            #a la fin d une sequence il est possible d avoir des GGG dans se cas
            # on verifie si la secore+1<len(line) car il ya un deuxieme G 
            #et 
            if (item < len(line) and (line[item]=="G" or line[item]=="g")):
                liste.append(1)
                #print liste
                if(item+1< len(line) and (line[item+1]=="G" or line[item+1]=="g")):
                    liste[item]=2
                    liste.append(2)
                    if (item+2< len(line) and (line[item+2]=="G" or line[item+2]=="g")):
                        liste[item+1]=3
                        liste[item]=3
                        liste.append(3)
                        if (item+3< len(line) and (line[item+3]=="G" or line[item+3]=="g")):
                            liste[item]=4
                            liste[item+1]=4
                            liste[item+2]=4
                            liste.append(4)
                            item=item+1
                        item=item+1
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="G" or line[item]=="g")):
                        liste.append(4)
                        item=item+1
        
            #elif (item < len(line) and (line[item]=="T" or line[item]=="A"  or line[item]=="t" or line[item]=="a" or line[item]=="U"or line[item]=="u" or
            #line[item]=="-" or line[item]=="N" or line[item]=="_" or line[item]=="Y"
            #                            or line[item]=="W" or line[item]=="R" or
            #                           line[item]=="K" or line[item]=="M"or #line[item]=="S" or line[item]=="B"
            #                       or line[item]=="V"or line[item]=="D"or line[item]=="H"or line[item]=="N")):
            #   liste.append(0)
            #   item=item+1
            elif (item < len(line) and line[item]!="G" and line[item]!="g" and line[item]!= "C" and line[item]!="c" ):
                        liste.append(0)
                        item=item+1
                
            elif(item < len(line) and (line[item]=="C" or line[item]=="c")):
                liste.append(-1)
                if(item+1< len(line) and (line[item+1]=="C" or line[item+1]=="c" )):
                    liste[item]=-2
                    liste.append(-2)
                    if (item+2< len(line) and (line[item+2]=="C" or line[item+2]=="c" )):
                        liste[item+1]=-3
                        liste[item]=-3
                        liste.append(-3)
                        if (item+3< len(line) and (line[item+3]=="C" or line[item+3]=="c"  )):
                            liste[item]=-4
                            liste[item+1]=-4
                            liste[item+2]=-4
                            liste.append(-4)
                            item=item+1
                        item=item+1   
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="C" or line[item]=="c")):
                    liste.append(-4)
                    item=item+1
            
            else:
                    item=item+1 #la fin du la ligne ou y a des entrers
        return line, liste
 
    def CalScore(self,liste, k):
        Score_Liste=[]
        #calcule de la moynne des scores pour toutes les sequences - les derniers k bases
        for i in range (len (liste)-(k-1)):
            #print (len(liste)), i, k
            j,Sum=0,0
            while (j<k):
                #print j, i
                Sum=Sum+liste[i]
                j=j+1
                i=i+1
            Mean=Sum/float(k)
            Score_Liste.append(Mean) 
        return Score_Liste
    
    ###################################################
    ###################################################
    def plot2(self,liste, repert, i):
        # make a little extra space between the subplots
        plt.subplots_adjust(wspace=1.0)
        dt = 1
        t = np.arange(0, len(liste), dt)
        figure= plt.figure()
        plt.plot(t, liste, 'b-')
        plt.xlim(0,len(liste))
        #figure.suptitle('Score of window nucleotides', fontsize=16)
        plt.xlabel('Position (ntS)')
        plt.ylabel('Score')
        plt.grid(True)
        figure.savefig(repert+'Score_plot_'+i+'.pdf', dpi=figure.dpi)
         
    """ 
    ###################################################
    ###################################################
    """               
    def GetG4(self,line,fileout,liste, Window,k, header, Len ):
        LG4=[]
        SEQ=">"+header+"\n Start \t End \t Sequence\t Length \t Score\n"
        fileout.write(SEQ)
        for i in range(len(liste)) :
            if (liste[i]>= float(Window) or liste[i]<= - float(Window)):
                seq=line[i:i+k]
                LG4.append(i)
                self.Write(fileout, i, k ,0,0, seq , k, liste[i])
                fileout.write("\n")
        return LG4
   
    def WriteSeq(self,line,fileout,liste, LISTE,header, F, Len ):
        i,k,I=0,0,0
        a=b=LISTE[i]
        MSCORE=[]
        SEQ=">"+header+"\nStart\tEnd\tSequence\tLength\tScore\tNBR\n"
        fileout.write(SEQ)
        if (len(LISTE)>1):
            c=LISTE[i+1]
            while (i< len(LISTE)-2):
                if(c==b+1):
                    k=k+1
                    i=i+1
                else:
                    I=I+1
                    seq=line[a:a+F+k]
                    sequence,liste2=self.BaseScore(seq)
                    self.Write(fileout, a, k ,F,0, seq ,len(seq) , round(np.mean(liste2),2))
                    MSCORE.append(abs(round(np.mean(liste2),2)))
                    fileout.write("\n")  
                    k=0
                    i=i+1
                    a=LISTE[i]
                b=LISTE[i] 
                c=LISTE[i+1] 
            I=I+1
            seq=line[a:a+F+k+1]
            sequence,liste2=self.BaseScore(seq)
            self.Write(fileout, a, k ,F,1, seq ,len(seq) , round(np.mean(liste2),2))
            MSCORE.append(abs(round(np.mean(liste2),2)))
            fileout.write("\t")
            fileout.write(str(I))
            fileout.write("\n") 
        #dans le cas ou il ya une seul sequence donc pas de chevauchement
        else:
            I=I+1
            seq=line[a:a+F]
            self.Write(fileout, a, 0 ,F,0, seq ,len(seq) , liste[a])
            MSCORE.append(abs(liste[a]))
            fileout.write("\t")
            fileout.write(str(I))
            fileout.write("\n")         
        return MSCORE
    
    def Write(self,fileout, i, k ,F,X, seq ,long, score):
        LINE=str(i)+" \t "+str(i+k+F+X)+" \t "+str(seq)+" \t "+str(long)+" \t "+str(score)
        fileout.write(LINE)
       
    #Len dans le cas ou la sequence fini avec des ----- donc il yaura une erreur



if __name__ == "__main__":
    try:
        inputfile, outputfile , window, score = main(sys.argv[1:])
        fname=inputfile.split("/")[-1]
        name=fname.split(".")
    #DIR="Results_"+str(name[0])
    #print DIR
    except ValueError:
        print '\033[1m' +"\n \t Oops! invalide parameters  \n" +'\033[0;0m'
        print "--------------------------------------------------------------------\n"
        sys.exit()
    except UnboundLocalError:
        print '\033[1m' +"\n \t Oops! invalide parameters  \n" +'\033[0;0m'
        print "--------------------------------------------------------------------\n"
        sys.exit()

    OPF= os.listdir(outputfile)

    flag=False
    if len(OPF)==0:
        flag=False
        DIR="Results_"+str(name[0])
    else:
        for dir in OPF:
            DIR="Results_"+str(name[0])
            if dir== DIR:
                flag=True
    if flag==True:
        shutil.rmtree(outputfile+"/"+DIR+"/")
        os.makedirs(outputfile+"/"+DIR+"/", mode=0777)        #
        print '\033[1m' +"\n \t Re-evaluation of G-quadruplex propensity with G4Hunter " +'\033[0;0m'
        print "\n#####################################"
        print "#    New Results directory Created  #"
        print "#####################################\n"
    else:
        os.makedirs(outputfile+"/"+DIR+"/", mode=0777)        #
        print "\n########################################################################"
        print "#                            Results directory Created                 #"
        print "########################################################################\n"
    
    #================================================================
    plot=[]
    fname=inputfile.split("/")[-1]
    filefasta=fname.split(".")
    filein=open(inputfile,"r")
    print "\n Input file:", '\033[1m' + filefasta[0]+'\033[0;0m'
    #repertoire des fichiers de sortie

    Res1file= open (outputfile +"/"+DIR+"/"+filefasta[0]+"-W"+ str(window)+"-S"+str(score)+".txt", "w")
    Res2file= open (outputfile +"/"+DIR+"/"+filefasta[0]+"-Merged.txt", "w")
    #=========================================
    
    startTime = time.time()
    
    
    soft1=Soft()
    ScoreListe, DNASeq, NumListe, HeaderListe=soft1.GFinder(filein, window)
    for i in range(len(DNASeq)):
        G4Seq=soft1.GetG4(DNASeq[i],Res1file, ScoreListe[i], float(score), int(window), HeaderListe[i], len(NumListe[i]))
        if (len(G4Seq)>0):
            MSCORE=soft1.WriteSeq(DNASeq[i],Res2file,ScoreListe[i], G4Seq, HeaderListe[i], int(window), len(NumListe[i]))
    #plot.append(MSCORE)

    malist, alllist=[], []
    for jj in range (len(ScoreListe[0])):
        cc, mean=0, 0
        for kk in range(len(ScoreListe)):
            cc+=ScoreListe[kk][jj]
        mean=cc/len(ScoreListe)
        alllist.append(mean)
        if abs(mean) >=score :
            malist.append(mean)
        else:
            malist.append(0)

    #soft1.plot2(ScoreListe[0], outputfile +"/Results/")
    #soft1.plot2(malist, outputfile +"/"+DIR+"/", "sc")
    soft1.plot2(alllist, outputfile +"/"+DIR+"/", "all")
    filein.close()
    fin=time.time()

    print "\n Results files and Score Figure are created in:   "#,fin-startTime, "secondes"
    print '\033[1m' + outputfile,DIR,"/","\n "+'\033[0;0m'


    Res1file.close()
    Res2file.close()

    
    
