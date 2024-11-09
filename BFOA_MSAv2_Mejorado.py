from bacteria import bacteria
from chemiotaxis import chemiotaxis
import csv
import os
import numpy

poblacion = []
path = 'multiFasta.fasta'
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 25
tumbo = 2                                            #numero de gaps a insertar 
nado = 1
chemio = chemiotaxis()
veryBest = bacteria(path)                #mejor bacteria   
tempBacteria = bacteria(path)            #bacteria temporal para validaciones
original = bacteria(path)                #bacteria original sin gaps
globalNFE = 0      #numero de evaluaciones de la funcion objetivo

dAttr= 1.2 #0.1
wAttr= 0.1 #0.2
hRep=dAttr
wRep= 10    #10


def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction
    
def validaSecuencias(path, veryBest):  
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-","")
   

    #valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return


for i in range(numeroDeBacterias):                                            #poblacion inicial
    poblacion.append(bacteria(path))

with open('resultados2.csv', mode='a', newline='') as file:
    writer = csv.writer(file)
    if os.path.getsize('resultados2.csv') == 0:
        writer.writerow(["Mejor Fitness", "NFE"])

    for _ in range(iteraciones):                                                  #numero de iteraciones  
        for bacteria in poblacion:
            bacteria.tumboNado(tumbo)
            bacteria.tumboNado(nado)
            bacteria.autoEvalua()  
        chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)                 #d_attr, w_attr, h_rep, w_rep):
        globalNFE += chemio.parcialNFE 
        best = max(poblacion, key=lambda x: x.fitness)
        if (veryBest == None) or (best.fitness > veryBest.fitness):
            clonaBest(veryBest, best)
        
        chemio.eliminarClonar(path, poblacion)
        chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)                #inserta  bacterias aleatorias


    writer.writerow([ veryBest.fitness, globalNFE])

veryBest.showGenome()
validaSecuencias(path, veryBest)
