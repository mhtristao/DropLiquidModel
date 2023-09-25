"""
@Mario Tristao, 2023.
Neutron Star Introduction (IFES)
------------------------------------------------------------------
This program simulate a semi-empirical mass formula represented by
liquid-drop-model proposed by George Gamow wich treat nucleus as a
drop of incompressible fluid of very high density.
M_(Z,A) = f_0 + f_1 + f_2 + f_3 + f_4 + f_5  
------------------------------------------------------------------
REF.:
    https://en.wikipedia.org/wiki/Semi-empirical_mass_formula#The_formula
    @book - Quantum Physics by Eisberg and Resnick, ch. 15.5
"""
import matplotlib.pyplot as plt
import numpy as np

#parametros
a1 = 0.01691
a2 = 0.01911
a3 = 0.000763
a4 = 0.10175
a5 = 0.012

#N = A - Z #Numero de Neutron
def massa(Z,A): #Massa atomica nuclear
    #1term. massa de um H_1 'u' 2term. Massa de um Neutron 'u'(UMA)
    #correcoes devido a energia de ligacao... 
    corF0 = 1.007825*Z + 1.008665*(A-Z) 
    volF1 =-a1 * A  #contribuicoes do volume
    supF2 = a2 * A**(2/3) #contribuicoes de superficie
    couF3 = a3 * Z**2/A**(1/3) #contribuicoes coulombianas
    assF4 = a4 * ((Z - (A/2))**2)/A #contribuicoes de assimetria 
    zParidade = Z % 2 #contribuição de emparelhamento - TERMOS DE PARIDADE
    nParidade = (A-Z) % 2
    if zParidade == 0 and nParidade == 0:
        empF5 = -a5 / np.sqrt(A)
    elif zParidade == 0 and nParidade !=0:
        empF5 = -a5 / np.sqrt(A)
    elif zParidade != 0 and nParidade !=0:
        empF5 = a5 / np.sqrt(A)
    elif zParidade !=0 and nParidade !=0:
        empF5 = a5 / np.sqrt(A)
    else:
        return 0
    return corF0 + volF1 + supF2 + couF3 + assF4 + empF5

#Energia de Ligacao - Transformação de UMA -> MeV
def enerL(Z,A):
    corF0 = 1.007825*Z + 1.008665*(A-Z)
    MeV1 = (corF0 - massa(Z,A)) * 931.5
    return MeV1

#Energia por Nucleon - MeV
def enerNu(Z,A):
    enerNucle = enerL(Z,A)/A
    return enerNucle

xModel = []
yModel = []
for i in range(250):
    i = i + 1
    z = i / 2.0
    xModel.append(i)
    yModel.append(enerNu(z,i))
    
data = np.array([xModel,yModel]).transpose()
np.savetxt('table1.dat', data, fmt=['%.1f\t', '%.4f'], header = 'xModel    yModel', comments='#')

model = plt
model.plot(xModel, yModel,'.r',label = 'Model')
model.plot(xModel, yModel,'-b')
model.title('Modelo Nuclear Gota-Líquida')
model.xlim([0,260])
model.ylim([0,10])
model.xlabel("Massa Nuclear A")
model.ylabel("Energia de Ligacao / nucleon (MeV)")
model.savefig("plot.png")
model.show()

print('Teste\n')
print('Cu Z=29 Massa Nuclear do Isotopo A=64')
print(massa(29,64))
print(enerL(29,64))
print(enerNu(29,64))