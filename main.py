from kivy.config import Config
Config.set('kivy', 'keyboard_mode', 'systemanddock')
Config.set('kivy', 'keyboard_layout', 'keypad.json')

import sys
sys.setrecursionlimit(50000)

import sympy
from chemhelper import compounds, conversions
from kivy.lang import Builder
from kivy.properties import ObjectProperty
from kivy.uix.screenmanager import Screen
from kivymd.app import MDApp
from kivymd.uix.boxlayout import BoxLayout
from math import pi
from trianglesolver import solve
from kivymd.uix.picker import MDThemePicker
from kivymd.uix.dialog import MDDialog
import re
from sympy import Matrix, lcm


# MORE FEATURES BEFORE INITIAL RELEASE:
#   CHEM: ***fix redox rxns
#       Calorimetry
#       Empirical/Molecular Formula
#       Percent Yield
#       Hydrates
#       Enthalpy/Entropy/Gibbs Free Energy
#       Oxidation Numbers***** NEED TO DO IONS
#       Solution Stoichiometry--chemhelper module TABS?


def decode(string="", addbefore=""):
	for i in string:
	    for j in superscripts.values():
	        if i == j:
	            replacement = addbefore + list(superscripts.keys())[list(superscripts.values()).index(j)]
	            string = string.replace(i, replacement)
	for i in string:
	    for j in subscripts.values():
	        if i == j:
	            replacement = addbefore + list(subscripts.keys())[list(subscripts.values()).index(j)]
	            string = string.replace(i, replacement)
	return string
def uncode(string="", which=""):
	for i in string:
	    if i.isdigit() and 'subscripts' in which:
	        string = string.replace(i, subscripts[i])
	    if (i.isdigit() or i == '+' or i == '-') and 'superscripts' in which and (i[i.index(i)-1] == '⁺' or i[i.index(i)-1] == '\u207B'):
	        string = string.replace(i, superscripts[i])
	return string


def bal(equation):
	elementList=[]
	elementMatrix=[]
	equation = equation.replace('+', ' ')
	r, p = equation.split('➞')
	reactants=r.split()
	products=p.split()
	
	def addToMatrix(element, index, count, side):
		if(index == len(elementMatrix)):
			elementMatrix.append([])
			for x in elementList:
				elementMatrix[index].append(0)
		if(element not in elementList):
			elementList.append(element)
			for i in range(len(elementMatrix)):
				elementMatrix[i].append(0)
		column=elementList.index(element)
		elementMatrix[index][column]+=count*side
    
	def findElements(segment,index, multiplier, side):
	    elementsAndNumbers=re.split('([A-Z][a-z]?)',segment)
	    i=0
	    while(i<len(elementsAndNumbers)-1):#last element always blank
	          i+=1
	          if(len(elementsAndNumbers[i])>0):
	            if(elementsAndNumbers[i+1].isdigit()):
	                count=int(elementsAndNumbers[i+1])*multiplier
	                addToMatrix(elementsAndNumbers[i], index, count, side)
	                i+=1
	            else:
	                addToMatrix(elementsAndNumbers[i], index, multiplier, side)        
    
	def compoundDecipher(compound, index, side):
	    segments=re.split('(\([A-Za-z0-9]*\)[0-9]*)',compound)    
	    for segment in segments:
	        if segment.startswith("("):
	            segment=re.split('\)([0-9]*)',segment)
	            multiplier=int(segment[1])
	            segment=segment[0][1:]
	        else:
	            multiplier=1
	        findElements(segment, index, multiplier, side)
	            
	for i in range(len(reactants)):
	    compoundDecipher(reactants[i],i,1)
	for i in range(len(products)):
	    compoundDecipher(products[i],i+len(reactants),-1)
	elementMatrix = Matrix(elementMatrix)
	elementMatrix = elementMatrix.transpose()
	solution=elementMatrix.nullspace()[0]
	multiple = lcm([val.q for val in solution])
	solution = multiple*solution
	coEffi=solution.tolist()
	output=""
	for i in range(len(reactants)):
	    output+=str(coEffi[i][0])+uncode(reactants[i])
	    if i<len(reactants)-1:
	       output+=" + "
	output+=" ➞ "
	for i in range(len(products)):
	   output+=str(coEffi[i+len(reactants)][0])+uncode(products[i])
	   if i<len(products)-1:
	       output+=" + "

	eq = {}
	for i in range(len(reactants)):
		eq[reactants[i]] = coEffi[i][0]
	for i in range(len(products)):
		eq[products[i]] = coEffi[i+len(reactants)][0]

	result = []
	result.append(output)
	result.append(eq)
	return result


symbols = {"hydrogen": 'H', "helium": 'He', "lithium": 'Li', 'beryllium': 'Be', 'boron': 'B',
       "carbon": 'C', "nitrogen": 'N', "oxygen": 'O', "fluorine": 'F', "neon": 'Ne',
       "sodium": 'Na', "magnesium": 'Mg', "aluminium": 'Al', "silicon": 'Si', "phosphorus": 'P',
       'sulfur': 'S', 'chlorine': 'Cl', "argon": 'Ar', 'potassium': 'K', 'calcium': 'Ca',
       'scandium': 'Sc', 'titanium': 'Ti', 'vanadium': 'V', 'chromium': 'Cr', 'manganese': 'Mn',
       'iron': 'Fe', 'cobalt': 'Co', 'nickel': 'Ni', 'copper': 'Cu', 'zinc': 'Zn', 'gallium': 'Ga',
       'germanium': 'Ge', 'arsenic': 'As', 'selenium': 'Se', 'bromine': 'Br', 'krypton': 'Kr',
       'rubidium': 'Rb', 'strontium': 'Sr', 'yttrium': 'Y', 'Zirconium': 'Zr', 'niobium': 'Nb',
       'molybdenum': 'Mo', 'technetium': 'Tc', 'ruthenium': 'Ru', 'rhodium': 'Rh', 'palladium': 'Pd',
       'silver': 'Ag', 'cadmium': 'Cd', 'indium': 'In', 'tin': 'Sn', 'antimony': 'Sb',
       'tellurium': 'Te', 'iodine': 'I', 'xenon': 'Xe', 'caesium': 'Cs', 'barium': 'Ba',
       'lanthanum': 'La', 'cerium': 'Ce', 'praseodymium': 'Pr', 'neodymium': 'Nd', 'promethium': 'Pm',
       'samarium': 'Sm', 'europium': 'Eu', 'gadolinium': 'Gd', 'terbium': 'Tb', 'dysprosium': 'Dy',
       'holmium': 'Ho', 'erbium': 'Er', 'thulium': 'Tm', 'ytterbium': 'Yb', 'lutetium': 'Lu',
       'hafnium': 'Hf', 'tantalum': 'Ta', 'tungsten': 'W', 'rhenium': 'Re', 'osmium': 'Os',
       'iridium': 'Ir', 'platinum': 'Pt', 'gold': 'Au', 'mercury': 'Hg', 'thallium': 'Tl',
       'lead': 'Pb', 'bismuth': 'Bi', 'polonium': 'Po', 'astatine': 'At', 'radon': 'Rn',
       'francium': 'Fr', 'radium': 'Ra', 'actinium': 'Ac', 'thorium': 'Th', 'protactinium': 'Pa',
       'uranium': 'U', 'neptunium': 'Np', 'plutonium': 'Pu', 'americium': 'Am', 'curium': 'Cm',
       'berkelium': 'Bk', 'californium': 'Cf', 'einsteinium': 'Es', 'fermium': 'Fm', 'mendelevium': 'Md',
       'nobelium': 'No', 'lawrencium': 'Lr', 'rutherfordium': 'Rf', 'dubnium': 'Db', 'seaborgium': 'Sg',
       'bohrium': 'Bh', 'hassium': 'Hs', 'Meitnerium': 'Mt', 'darmstadtium': 'Ds', 'roentgenium': 'Rg',
       'copernicium': 'Cn', 'nihonium': 'Nh', 'flerovium': 'Fl', 'moscovium': 'Mc', 'livermorium': 'Lv',
       'tennessine': 'Ts', 'Oganesson': 'Og'
       }

ele_names = list(symbols.keys())
symbols_list = list(symbols.values())

symbol_valences = {'H':(1,-1),'Li':(1),'Zn':(2),'Sb':(3,5),'Bi':(3,5),'Ag':(1),'Au':(1,3),
        'Sn':(2,4),'Hg':(1,2),'Ni':(2,3),'Pb':(2,4),'Cu':(1,2),
        'Co':(2,3),'Mn':(2,3),'Cr':(2,3),'Na':(1),'K':(1),'Rb':(1),
        'Cs':(1),'Fr':(1),'Be':(2),'Mg':(2),'Ca':(2),'Sr':(2),'Ba':(2),
        'Ra':(2),'B':(3),'Al':(3),'Ga':(3),'In':(3),'Tl':(3),'C':(4,-4),
        'Si':(4,-4),'Ge':(4,-4),'N':(-3),'P':(-3),'As':(-3),'O':(-2),
        'S':(-2),'Se':(-2),'Te':(-2),'Po':(-2),'F':(-1),'Cl':(-1),'Br':(-1),
        'I':(-1),'At':(-1)
        } # add more

polyatomics = {'NH4':+1,'H3O':+1,'C2H3O2':-1,
           'CN':-1,'H2PO4':-1,'HCO3':-1,
           'HSO4':-1,'OH':-1,'ClO':-1,'ClO2':-1,
           'ClO3':-1,'ClO4':-1,'BrO':-1,'BrO2':-1,'BrO3':-1,
           'BrO4':-1,'IO':-1,'IO2':-1,'IO3':-1,'IO4':-1,
           'NO3':-1,'NO2':-1,'MnO4':-1,'SCN':-1,
           'CO3':-2,'CrO4':-1,'Cr2O7':-2,'HPO4':-2,
           'SO4':-2,'SO3':-2,'S2O3':-2,'PO4':-3,'PO3':-3,'AsO4':-3
          }


polyatomics_keys = list(polyatomics.keys())
polyatomics_values = list(polyatomics.values())

# print(ele_names[symbols_list.index('Cn')])

sqrt_symbol = u'\u221A'
one_symbol = u'\u00B9'
squared_symbol = u'\u00B2'
cubed_symbol = u'\u00B3'
four_symbol = u'\u2074'
five_symbol = u'\u2075'
six_symbol = u'\u2076'
seven_symbol = u'\u2077'
eight_symbol = u'\u2078'
nine_symbol = u'\u2079'
ten_symbol = one_symbol + u'\u2070'
eleven_symbol = one_symbol + one_symbol
twelve_symbol = one_symbol + squared_symbol
thirteen_symbol = one_symbol + cubed_symbol
fourteen_symbol = one_symbol + four_symbol

subscripts = {'1':u'\u2081','2':u'\u2082','3':u'\u2083','4':u'\u2084','5':u'\u2085',
        '6':u'\u2086','7':u'\u2087','8':u'\u2088','9':u'\u2089','0':u'\u2080'}

superscripts = {'1':one_symbol,'2':squared_symbol,'3':cubed_symbol,'4':four_symbol,'5':five_symbol,
        '6':six_symbol,'7':seven_symbol,'8':eight_symbol,'9':nine_symbol,'0':u'\u2070','+':'⁺','-':'\u207B'}

def get_element_symbol(self, name):
	return symbols_list[ele_names.index(name.lower())]

elec_configs = {'H':('1s'+one_symbol+'','1s'+one_symbol+''),'He':('1s'+squared_symbol+'','1s'+squared_symbol+''),'Li':('1s'+squared_symbol+'2s'+one_symbol+'','[He]2s'+one_symbol+''),'Be':('1s'+squared_symbol+'2s'+squared_symbol+'','[He]2s'+squared_symbol+''),
            'B':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+one_symbol+'','[He]2s'+squared_symbol+'2p'+one_symbol+''),'C':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+squared_symbol+'','[He]2s'+squared_symbol+'2p'+squared_symbol+''),'N':('1s'+squared_symbol+'2s'+squared_symbol+'2d'+cubed_symbol+'','[He]2s'+squared_symbol+'2d'+cubed_symbol+''),
            'O':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+four_symbol+'','[He]2s'+squared_symbol+'2p'+four_symbol+''),'F':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+five_symbol+'','[He]2s'+squared_symbol+'2p'+five_symbol+''),'Ne':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'','[He]2s'+squared_symbol+'2p'+six_symbol+''),
            'Na':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+one_symbol+'','[Ne]3s'+one_symbol+''),'Mg':('1s'+one_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'','[Ne]3s'+squared_symbol+''),'Al':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+one_symbol+'','[Ne]3s'+squared_symbol+'3p'+one_symbol+''),
            'Si':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+squared_symbol+'','[Ne]3s'+squared_symbol+'3p'+squared_symbol+''),'P':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3d'+cubed_symbol+'','[Ne]3s'+squared_symbol+'3d'+cubed_symbol+''),'S':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+four_symbol+'','[Ne]3s'+squared_symbol+'3p'+four_symbol+''),
            'Cl':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+five_symbol+'','[Ne]3s'+squared_symbol+'3p'+five_symbol+''),'Ar':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'','[Ne]3s'+squared_symbol+'3p'+six_symbol+''),'K':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'4s'+one_symbol+'','[Ar]4s'+one_symbol+''),
            'Ca':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'4s'+squared_symbol+'','[Ar]4s'+squared_symbol+''),'Sc':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+one_symbol+'4s'+squared_symbol+'','[Ar]3d'+one_symbol+'4s'+squared_symbol+''),'Ti':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+squared_symbol+'4s'+squared_symbol+'','[Ar]3d'+squared_symbol+'4s'+squared_symbol+''),
            'V':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+cubed_symbol+'4s'+squared_symbol+'','[Ar]3d'+cubed_symbol+'4s'+squared_symbol+''),'Cr':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+five_symbol+'4s'+one_symbol+'','[Ar]3d'+five_symbol+'4s'+one_symbol+''),'Mn':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+five_symbol+'4s'+squared_symbol+'','[Ar]3d'+five_symbol+'4s'+squared_symbol+''),
            'Fe':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+six_symbol+'4s'+squared_symbol+'','[Ar]3d'+six_symbol+'4s'+squared_symbol+''),'Co':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+seven_symbol+'4s'+squared_symbol+'','[Ar]3d'+seven_symbol+'4s'+squared_symbol+''),'Ni':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+eight_symbol+'4s'+squared_symbol+'','[Ar]3d'+eight_symbol+'4s'+squared_symbol+''),
            'Cu':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+one_symbol+'','[Ar]3d'+ten_symbol+'4s'+one_symbol+''),'Zn':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'','[Ar]3d'+ten_symbol+'4s'+squared_symbol+''),'Ga':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+one_symbol+'','[Ar]3d'+ten_symbol+'4s'+squared_symbol+'4p'+one_symbol+''),
            'Ge':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+squared_symbol+'','[Ar]3d'+ten_symbol+'4s'+squared_symbol+'4p'+squared_symbol+''),'As':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4d'+cubed_symbol+'','[Ar]3d'+ten_symbol+'4s'+squared_symbol+'4d'+cubed_symbol+''),'Se':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+four_symbol+'','[Ar]3d'+ten_symbol+'4s'+squared_symbol+'4p'+four_symbol+''),
            'Br':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+five_symbol+'','[Ar]3d'+ten_symbol+'4s'+squared_symbol+'4p'+five_symbol+''),'Kr':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'','[Ar]3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+''),'Rb':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'5s'+one_symbol+'','[Kr]5s'+one_symbol+''),
            'Sr':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'5s'+squared_symbol+'','[Kr]5s'+squared_symbol+''),'Y':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+one_symbol+'5s'+squared_symbol+'','[Kr]4d'+one_symbol+'5s'+squared_symbol+''),'Zr':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+squared_symbol+'5s'+squared_symbol+'','[Kr]4d'+squared_symbol+'5s'+squared_symbol+''),
            'Nb':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+four_symbol+'5s'+one_symbol+'','[Kr]4d'+four_symbol+'5s'+one_symbol+''),'Mo':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+five_symbol+'5s'+one_symbol+'','[Kr]4d'+five_symbol+'5s'+one_symbol+''),'Tc':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+five_symbol+'5s'+squared_symbol+'','[Kr]4d'+five_symbol+'5s'+squared_symbol+''),
            'Ru':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+seven_symbol+'5s'+one_symbol+'','[Kr]4d'+seven_symbol+'5s'+one_symbol+''),'Rh':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+eight_symbol+'5s'+one_symbol+'','[Kr]4d'+eight_symbol+'5s'+one_symbol+''),'Pd':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'','[Kr]4d'+ten_symbol+''),
            'Ag':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+one_symbol+'','[Kr]4d'+ten_symbol+'5s'+one_symbol+''),'Cd':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'','[Kr]4d'+ten_symbol+'5s'+squared_symbol+''),'In':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+one_symbol+'','[Kr]4d'+ten_symbol+'5s'+squared_symbol+'5p'+one_symbol+''),
            'Sn':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+squared_symbol+'','[Kr]4d'+ten_symbol+'5s'+squared_symbol+'5p'+squared_symbol+''),'Sb':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5d'+cubed_symbol+'','[Kr]4d'+ten_symbol+'5s'+squared_symbol+'5d'+cubed_symbol+''),'T':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+four_symbol+'','[Kr]4d'+ten_symbol+'5s'+squared_symbol+'5p'+four_symbol+''),
            'I':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+five_symbol+'','[Kr]4d'+ten_symbol+'5s'+squared_symbol+'5p'+five_symbol+''),'Xe':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+six_symbol+'','[Kr]4d'+ten_symbol+'5s'+squared_symbol+'5p'+six_symbol+''),'Cs':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+one_symbol+'','[Xe]6s'+one_symbol+''),
            'Ba':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]6s'+squared_symbol+''),'La':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+one_symbol+'6s'+squared_symbol+'','[Xe]5d'+one_symbol+'6s'+squared_symbol+''),'Ce':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+one_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+one_symbol+'6s'+squared_symbol+'','[Xe]4f'+one_symbol+'5d'+one_symbol+'6s'+squared_symbol+''),
            'Pr':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4d'+cubed_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4d'+cubed_symbol+'6s'+squared_symbol+''),'Nd':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+four_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+four_symbol+'6s'+squared_symbol+''),'Pm':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+five_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+five_symbol+'6s'+squared_symbol+''),
            'Sm':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+six_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4d'+cubed_symbol+'6s'+squared_symbol+''),'Eu':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+seven_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+four_symbol+'6s'+squared_symbol+''),'Gd':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+seven_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+one_symbol+'6s'+squared_symbol+'','[Xe]4f'+seven_symbol+'5d'+one_symbol+'6s'+squared_symbol+''),
            'Tb':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+nine_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+nine_symbol+'6s'+squared_symbol+''),'Dy':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+ten_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+ten_symbol+'6s'+squared_symbol+''),'Ho':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+eleven_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+eleven_symbol+'6s'+squared_symbol+''),
            'Er':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+twelve_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+twelve_symbol+'6s'+squared_symbol+''),'Tm':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+thirteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+thirteen_symbol+'6s'+squared_symbol+''),'Yb':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'6s'+squared_symbol+''),
            'Lu':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+one_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+one_symbol+'6s'+squared_symbol+''),'Hf':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+squared_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+squared_symbol+'6s'+squared_symbol+''),'Ta':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+cubed_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+cubed_symbol+'6s'+squared_symbol+''),
            'W':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+four_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+four_symbol+'6s'+squared_symbol+''),'Re':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+five_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+five_symbol+'6s'+squared_symbol+''),'Os':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+six_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+six_symbol+'6s'+squared_symbol+''),
            'Ir':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+seven_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+seven_symbol+'6s'+squared_symbol+''),'Pt':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+nine_symbol+'6s'+one_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+nine_symbol+'6s'+one_symbol+''),'Au':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+one_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+one_symbol+''),
            'Hg':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+squared_symbol+''),'Tl':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+one_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+one_symbol+''),'Pb':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+squared_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+squared_symbol+''),
            'Bi':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6d'+cubed_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6d'+cubed_symbol+''),'Po':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+four_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+four_symbol+''),'At':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+five_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+five_symbol+''),
            'Rn':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+six_symbol+'','[Xe]4f'+fourteen_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+six_symbol+''),'Fr':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+one_symbol+'','[Rn]7s'+one_symbol+''),'Ra':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]7s'+squared_symbol+''),
            'Ac':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+one_symbol+'7s'+squared_symbol+'','[Rn]6d'+one_symbol+'7s'+squared_symbol+''),'Th':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+squared_symbol+'7s'+squared_symbol+'','[Rn]6d'+squared_symbol+'7s'+squared_symbol+''),'Pa':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+squared_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+one_symbol+'7s'+squared_symbol+'','[Rn]5f'+squared_symbol+'6d'+one_symbol+'7s'+squared_symbol+''),
            'U':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5d'+cubed_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+one_symbol+'7s'+squared_symbol+'','[Rn]5d'+cubed_symbol+'6d'+one_symbol+'7s'+squared_symbol+''),'Np':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+four_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+one_symbol+'7s'+squared_symbol+'','[Rn]5f'+four_symbol+'6d'+one_symbol+'7s'+squared_symbol+''),'Pu':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+six_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+six_symbol+'7s'+squared_symbol+''),
            'Am':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+seven_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+seven_symbol+'7s'+squared_symbol+''),'Cm':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+seven_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+one_symbol+'7s'+squared_symbol+'','[Rn]5f'+seven_symbol+'6d'+one_symbol+'7s'+squared_symbol+''),'Bk':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+nine_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+nine_symbol+'7s'+squared_symbol+''),
            'Cf':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+ten_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+ten_symbol+'7s'+squared_symbol+''),'Es':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+eleven_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+eleven_symbol+'7s'+squared_symbol+''),'Fm':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+twelve_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+twelve_symbol+'7s'+squared_symbol+''),
            'Md':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+thirteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+thirteen_symbol+'7s'+squared_symbol+''),'No':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'7s'+squared_symbol+''),'Lr':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'7s'+squared_symbol+'7p'+one_symbol+'','[Rn]5f'+fourteen_symbol+'7s'+squared_symbol+'7p'+one_symbol+''),
            'Rf':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+squared_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+squared_symbol+'7s'+squared_symbol+''),'Db':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+cubed_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+cubed_symbol+'7s'+squared_symbol+''),'Sg':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+four_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+four_symbol+'7s'+squared_symbol+''),
            'Bh':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+five_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+five_symbol+'7s'+squared_symbol+''),'Hs':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+six_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+six_symbol+'7s'+squared_symbol+''),'Mt':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+seven_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+seven_symbol+'7s'+squared_symbol+''),
            'Ds':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+eight_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+eight_symbol+'7s'+squared_symbol+''),'Rg':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+nine_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+nine_symbol+'7s'+squared_symbol+''),'Cn':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+ten_symbol+'7s'+squared_symbol+''),
            'Nh':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+one_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+one_symbol+''),'Fl':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+squared_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+squared_symbol+''),'Mc':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7d'+cubed_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7d'+cubed_symbol+''),
            'Lv':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+four_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+four_symbol+''),'Ts':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+five_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+five_symbol+''),'Og':('1s'+squared_symbol+'2s'+squared_symbol+'2p'+six_symbol+'3s'+squared_symbol+'3p'+six_symbol+'3d'+ten_symbol+'4s'+squared_symbol+'4p'+six_symbol+'4d'+ten_symbol+'4f'+fourteen_symbol+'5s'+squared_symbol+'5p'+six_symbol+'5d'+ten_symbol+'5f'+fourteen_symbol+'6s'+squared_symbol+'6p'+six_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+six_symbol+'','[Rn]5f'+fourteen_symbol+'6d'+ten_symbol+'7s'+squared_symbol+'7p'+six_symbol+'')
}


class ContentNavigationDrawer(BoxLayout):
	screen_manager = ObjectProperty()
	nav_drawer = ObjectProperty()
	search = ObjectProperty()

	def text_search(self):
	    search_text = self.search.text
	    self.search.text = ""
	    for i in self.ids.keys():
	        if search_text not in i:
	            self.ids.mylist.remove_widget(self.ids[i])
	    # Trouble getting MDList to come back after search
	    # for i in self.ids.keys():
	    #     self.ids.mylist.add_widget(self.ids[i])

class ContentNavigationDrawer(BoxLayout):
	screen_manager = ObjectProperty()
	nav_drawer = ObjectProperty()
	search = ObjectProperty()

	def text_search(self):
	    search_text = self.search.text
	    self.search.text = ""
	    for i in self.ids.keys():
	        if search_text not in i:
	            self.ids.mylist.remove_widget(self.ids[i])
	    # Trouble getting MDList to come back after search
	    # for i in self.ids.keys():
	    #     self.ids.mylist.add_widget(self.ids[i])

class HomeScreen(Screen):
	pass
class Chemistry(Screen):
	pass
class ShorthandNotation(Screen):
	pass
class MolarMass(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def mass(self):
	    try:
	        self.data.text = decode(self.data.text)
	        self.result.text = str(compounds.mass(self.data.text)) + ' grams/mole'
	        self.data.text = ""
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""
class PercentComposition(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def composition(self):
	    try:
	        self.result.text = "Percents:\n" + str(compounds.percent_composition(decode(self.data.text))).replace(', dtype: float64', '')
	        self.data.text = ""
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""
class Molarity(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def calc(self):
	    try:
	        data_text = self.data.text
	        self.data.text = ""
	        data_list = data_text.split(',')

	        element = ''
	        m_data_text = data_text
	        if 'chemical' in m_data_text:
	            for i in m_data_text.split(','):
	                if 'chemical' in i:
	                    element = decode(i[i.index('=')+1:])

	        newstr = ''.join((ch if ch in '0123456789.' else ' ') for ch in m_data_text)
	        vals = [float(i) for i in newstr.split()]

	        elg_args = ['volume', 'mass', 'chemical', 'moles']
	        arguments = [None, None, None, None]
	        k = 0
	        for i in data_list:
	            for j in elg_args:
	                if j in i:
	                    if j == 'chemical':
	                        arguments.insert(elg_args.index(j), element)
	                    else:
	                        arguments.insert(elg_args.index(j), vals[k])
	                    k += 1

	        if data_list[0][data_list[0].index('=')+len(str(arguments[0]))+1:] == 'mL':
	            arguments[0] = arguments[0] / 1000

	        self.result.text = str(conversions.molarity(volume=arguments[0], mass=arguments[1], chemical=decode(arguments[2]), moles=arguments[3])) + ' M'
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""
class IdealGasLaw(Screen):
	data = ObjectProperty(None)
	volume_units = ObjectProperty(None)
	temp_units = ObjectProperty(None)
	pressure_units = ObjectProperty(None)
	result = ObjectProperty(None)

	item_index_a = 0
	item_index_b = 0
	item_index_c = 0

	vol_units = ['L','mL']
	temper_units = ['K','C','F']
	press_units = ['atm','kPa','mmHg','torr']


	def change_item_a(self):
	    index = self.item_index_a
	    self.item_index_a += 1
	    if index > len(self.vol_units) - 1:
	        index = 0
	        self.item_index_a = 1
	    self.volume_units.text = self.vol_units[index]

	def change_item_b(self):
	    index = self.item_index_b
	    self.item_index_b += 1
	    if index > len(self.temper_units) - 1:
	        index = 0
	        self.item_index_b = 1
	    self.temp_units.text = self.temper_units[index]

	def change_item_c(self):
	    index = self.item_index_c
	    self.item_index_c += 1
	    if index > len(self.press_units) - 1:
	        index = 0
	        self.item_index_c = 1
	    self.pressure_units.text = self.press_units[index]

	# haven't started method
	def ideal_gas(self):
	    try:
	        data_text = self.data.text
	        self.data.text = ""
	        vunits = self.volume_units.text
	        tunits = self.temp_units.text
	        punits = self.pressure_units.text
	        self.volume_units.text = 'Volume Units'
	        self.temp_units.text = 'Temperature Units'
	        self.pressure_units.text = 'Pressure Units'
	        self.item_index_a = 0
	        self.item_index_b = 0
	        self.item_index_c = 0

	        if '°' in data_text:
	            data_text = data_text.replace('°', '')

	        element = ''
	        if 'chemical' in data_text:
	            for i in data_text.split(','):
	                if 'chemical' in i:
	                    molar_m = compounds.mass(decode(i[i.index('=')+1:]))
	                    data_text = data_text.replace(i, '')

	        possible = ['volume','moles','temperature','pressure','mass','molar_mass','substance','solve_for']
	        entries = {}
	        for i in data_text.split(','):
	            for j in possible:
	                if j in i and (j != 'solve_for' and j != 'substance'):
	                    entries[j] = float(i[i.index('=')+1:])
	                elif j in i and (j == 'solve_for' or j == 'substance'):
	                    entries[j] = i[i.index('=')+1:]

	        solve = entries['solve_for'].lower()

	        # Volume Conversions
	        if solve != 'v' and solve not in 'volume' and vunits == 'mL':
	            entries['volume'] = entries['volume'] / 1000

	        # Temperature Conversions
	        if solve != 't' and solve not in 'temperature' and tunits == 'C':
	            entries['temperature'] += 273
	        elif solve != 't' and solve not in 'temperature' and tunits == 'F':
	            entries['temperature'] = (entries['temperature'] - 32) * (5/9) + 273

	        R = 0.0821
	        # R Determination
	        if punits == 'kPa':
	            R = 8.31
	        elif punits == 'mmHg' or punits == 'torr':
	            R = 62.4

	        if solve == 'p':
	            if 'mass' not in entries.keys():
	                self.result.text = str(entries['moles'] * R * entries['temperature'] / entries['volume']) + ' ' + punits
	                return
	            else:
	                if 'molar_mass' in entries.keys():
	                    self.result.text = str((entries['mass']/entries['molar_mass'])*R*entries['temperature']/entries['volume']) + ' ' + punits
	                    return
	                else:
	                    self.result.text = str((entries['mass']/molar_m)*R*entries['temperature']/entries['volume']) + ' ' + punits
	                    return
	        elif solve == 'v':
	            if 'mass' not in entries.keys():
	                self.result.text = str(entries['moles']*R*entries['temperature']/entries['pressure']) + ' L'
	            else:
	                if 'molar_mass' in entries.keys():
	                    self.result.text = str((entries['mass']/entries['molar_mass'])*R*entries['temperature']/entries['pressure']) + ' L'
	                else:
	                    self.result.text = str((entries['mass']/molar_m)*R*entries['temperature']/entries['pressure']) + ' L'
	        elif solve == 'n':
	            self.result.text = str((entries['pressure']*entries['volume'])/(R*entries['temperature'])) + ' moles'
	            return
	        elif solve == 't':
	            if 'mass' not in entries.keys():
	                self.result.text = str((entries['pressure']*entries['volume'])/(entries['moles']*R)) + ' K'
	            else:
	                if 'molar_mass' in entries.keys():
	                    self.result.text = str((entries['pressure']*entries['volume'])/(entries['mass']/entries['molar_mass']*R)) + ' K'
	                else:
	                    self.result.text = str((entries['pressure'] * entries['volume'])/(entries['mass'] / molar_m * R)) + ' K'
	        elif solve == 'r':
	            if punits == 'atm':
	                self.result.text = '0.0821 atm*L/mol*K'
	                return
	            elif punits == 'kPa':
	                self.result.text = '8.31 kPa*L/mol*K'
	                return
	            elif punits == 'mmHg':
	                self.result.text = '62.4 mmHg*L/mol*K'
	                return
	            elif punits == 'torr':
	                self.result.text = '62.4 torr*L/mol*K'
	                return
	        elif solve == 'm':
	            if 'molar_mass' in entries.keys():
	                self.result.text = str((entries['molar_mass']*entries['pressure']*entries['volume'])/(R*entries['temperature'])) + ' grams'
	                return
	            else:
	                self.result.text = str((molar_m * entries['pressure'] * entries['volume'])/(R * entries['temperature'])) + ' grams'
	                return
	        elif solve == 'mm':
	            self.result.text = str((entries['mass']*R*entries['temperature'])/(entries['pressure']*entries['volume'])) + ' grams/mole'
	            return

	        if vunits == 'mL':
	            self.result.text = str(float(self.result.text.split()[0])*1000) + ' mL'
	            return
	        if tunits == 'C':
	            self.result.text = str(float(self.result.text.split()[0])-273) + ' °C'
	            return
	        if tunits == 'F':
	            self.result.text = str((float(self.result.text.split()[0])-273)*(9/5)+32) + ' °F'
	            return
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""
	        self.volume_units.text = "Volume Units"
	        self.pressure_units.text = "Pressure Units"
	        self.temp_units.text = "Temperature Units"
class Solubility(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def soluble(self):
	    try:
	        data_text = decode(self.data.text)
	        self.data.text = ""
	        alkaline_metals = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
	        for i in alkaline_metals:
	            if i in data_text:
	                self.result.text = 'Soluble'
	                return
	        if 'NH4' in data_text:
	            self.result.text = 'Soluble'
	            return
	        if 'NO3' in data_text or 'ClO4' in data_text or 'ClO3' in data_text or 'C2H3O2' in data_text or 'CH3COO' in data_text:
	            self.result.text = 'Soluble'
	            return
	        if ('Cl' in data_text or 'Br' in data_text or 'I' in data_text) and ('Ag' not in data_text and 'Hg' not in data_text and 'Pb' not in data_text):
	            self.result.text = 'Soluble'
	            return
	        if 'SO4' in data_text and ('Hg' not in data_text and 'Pb' not in data_text and 'Sr' not in data_text and 'Ca' not in data_text and 'Ba' not in data_text):
	            self.result.text = 'Soluble'
	            return
	        if 'OH' in data_text and ('Ca' not in data_text and 'Sr' not in data_text and 'Ba' not in data_text):
	            self.result.text = 'Insoluble'
	            return
	        elif 'OH' in data_text and ('Ca' in data_text or 'Sr' in data_text or 'Ba' in data_text):
	            self.result.text = 'Soluble'
	            return
	        if 'PO4' in data_text or 'S' in data_text or 'CO3' in data_text or 'SO3' in data_text:
	            self.result.text = 'Insoluble'
	            return

	        self.result.text = 'Insoluble'
	        return
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""

class EmpMolFormula(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

# def calculate(self):
#     data_list = decode(self.data.text).split(',')
#     self.data.text = ""
#     entries = {}
#
#     possible = ['solve_for',]
#     num_ele = 0
#     syms_list = []
#     for i in data_list:
#         if 'num_elements' in i:
#             num_ele = int(i[i.index('=')+1:])
#
#     syms_list.append('' for i in range(num_ele))
#
#     index = 0
#     for i in data_list:
#         if 'percent' in i:
#             separate = i[i.index('=')+1:]
#             for j in separate:
#                 if j.isalpha():
#                     syms_list[index] += j
#             data_list[data_list.index(i)] = data_list[data_list.index(i)].replace(syms_list[index], '')
#             index += 1
#
#
#     entries = {}
#     k = 0
#     for i in range(num_ele):
#         possible.append('percent_'+str(i+1))
#         entries['element_'+str(i+1)] = syms_list[k]
#         k += 1
#
#
#     for i in data_list:
#         for j in possible:
#             if j in i and ('solve_for' in j):
#                 entries[j] = (i[i.index('=')+1:])
#             elif j in i and ('percent' in j):
#                 entries[j] = float(i[i.index('=')+1:])
#
#     ans = []
#     calcs = []
#     calcs.append(None for i in range(num_ele))
#     ans.append(None for i in range(num_ele))
#     if entries['solve_for'] == 'E':
#         for i in range(num_ele):
#             calcs[i] = round(entries['percent_'+str(i+1)] / compounds.mass(entries['element_'+str(i+1)]), 2)
#
#         test = False
#         for i in range(num_ele):
#             if str(calcs[i]/min(calcs))[2] == '9' or str(calcs[i]/min(calcs))[2] == '0':
#                 test = True
#             ans[i] = calcs[i] / min(calcs)
#
#         k = 2
#         bruh = True
#         if test:
#             while bruh:
#                 for i in range(num_ele):
#                     ans[i] *= 2
#                     if ans[i].is_integer():
#                         ans[i] = int(ans[i])
#                 k += 1



class Molality(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def calculate(self):
	    try:
	        data_text = decode(self.data.text)
	        self.data.text = ""
	        data_text = data_text.replace(' ', '')
	        desired = ['moles_solute','kgrams_solvent']
	        possible = ['solute','solvent','moles_solute','grams_solute','kgrams_solvent','grams_solvent','moles_solvent']
	        units = []
	        values = []
	        calculation = [None, None]
	        for i in data_text.split(','):
	            for j in possible:
	                if j == i[0:i.index('=')]:
	                    units.append(j)
	                    if j == 'solute' or j == 'solvent':
	                        values.append(compounds.mass(decode(i[i.index('=')+1:])))
	                    else:
	                        values.append(float(i[i.index('=')+1:]))

	        while desired != None:
	            for i in units:
	                index = units.index(i)
	                if i == 'moles_solute':
	                    if desired != None and 'moles_solute' in desired:
	                        desired = desired.remove('moles_solute')
	                    calculation[0] = values[index]
	                if i == 'kgrams_solvent':
	                    if desired != None and 'kgrams_solvent' in desired:
	                        desired = desired.remove('kgrams_solvent')
	                    calculation[1] = values[index]
	                if i == 'grams_solute':
	                    if desired != None and 'moles_solute' in desired:
	                        desired = desired.remove('moles_solute')
	                    molar_mass = values[units.index('solute')]
	                    calculation[0] = values[index] / molar_mass
	                if i == 'grams_solvent':
	                    if desired != None and 'kgrams_solvent' in desired:
	                        desired = desired.remove('kgrams_solvent')
	                    calculation[1] = values[index] / 1000
	                if i == 'moles_solvent':
	                    if desired != None and 'kgrams_solvent' in desired:
	                        desired = desired.remove('kgrams_solvent')
	                    molar_mass = values[units.index('solvent')]
	                    calculation[1] = values[index] * molar_mass / 1000

	        self.result.text = str(calculation[0] / calculation[1]) + ' molal'
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""
class ElectronConfigurations(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def configure(self):
	    try:
	        data_text = self.data.text.split(',')
	        self.data.text = ""
	        test = false
	        for i in symbols_list:
	            if data_text[0] == i:
	                test = true
	        if not test:
	            data_text[0] = get_element_symbol(self, data_text[0])

	        if 'true' in data_text[1]:
	            self.result.text = elec_configs[data_text[0]][1]
	        else:
	            self.result.text = elec_configs[data_text[0]][0]
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""

class OxidationNumbers(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def oxidate(self, extra, go):
	    try:
	        if not extra == None:
	            self.data.text = extra

	        self.result.text = ""
	        data_text = decode(self.data.text)
	        self.data.text = ""
	        has_charge = False
	        if '-' in data_text or '+' in data_text:
	            has_charge = True
	        ele_symbols_str = ''.join(filter(str.isalpha, data_text))
	        ele_symbols_list = []
	        for i in range(len(ele_symbols_str)):
	            if ele_symbols_str[i].isupper() and i != len(ele_symbols_str) - 1 and ele_symbols_str[i+1].islower():
	                ele_symbols_list.append(ele_symbols_str[i] + ele_symbols_str[i+1])
	            elif ele_symbols_str[i].isupper():
	                ele_symbols_list.append(ele_symbols_str[i])

	        poly_charge = 0
	        separate = ''
	        contains_poly = False
	        if go:
	            for i in polyatomics_keys:
	                if i in data_text:
	                    poly_charge = polyatomics[i]
	                    data_text = data_text.replace(i, '')
	                    contains_poly = True
	                    if poly_charge > 0:
	                        full = i + '+' + str(poly_charge)
	                        separate = self.oxidate(full, False)
	                    else:
	                        full = i+str(poly_charge)
	                        separate = self.oxidate(full, False)
	                    for j in symbols_list:
	                        for k in ele_symbols_list:
	                            if j == k and j in i:
	                                ele_symbols_list.remove(j)

	        digits = []
	        ox_nums = []
	        sign = ''
	        for i in range(len(ele_symbols_list)+1):
	            digits.append(None)
	        for i in data_text:
	            for j in ele_symbols_list:
	                if i.isdigit() and data_text.index(j) + len(j) == data_text.index(i):
	                    if data_text.index(i) != 0 and data_text[data_text.index(i)-1] != '+' and data_text[data_text.index(i)-1] != '-':
	                        # print(data_text[data_text.index(i) - 1])
	                        if digits[ele_symbols_list.index(j)] == None and None in digits:
	                            # print(data_text[data_text.index(i) - 1])
	                            digits.remove(None)
	                            digits.insert(ele_symbols_list.index(j),float(i))
	                            data_text = data_text.replace(i, '',1)
	                            break
	                if i.isdigit() and contains_poly:
	                    if None in digits:
	                        digits.remove(None)
	                    digits.append(float(i))
	                if i == '+' or i == '-':
	                    if i == '+':
	                        sign = '+'
	                    else:
	                        sign = '-'


	        charge = 0
	        if has_charge and data_text[data_text.index(sign)+1:] != '':
	            charge = int(data_text[data_text.index(sign)+1:])
	        elif has_charge:
	            charge = 1

	        if has_charge and len(ele_symbols_list) == 1:
	            self.result.text = ele_symbols_list[0] + ': ' + sign + str(charge)
	            return
	        if len(ele_symbols_list) == 1 and not contains_poly:
	            self.result.text = ele_symbols_list[0] + ': 0'
	            return

	        alkaline_metals = ['Li','Na','K','Rb','Cs','Fr']
	        alkaline_earth_metals = ['Be','Mg','Ca','Sr','Ba','Ra']

	        def contains_metal(self):
	            metals = alkaline_metals + alkaline_earth_metals
	            metals += ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
	                       'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
	                       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',
	                       'Lr','Rf','Db','Sg','Bh','Hs','Mt','La','Ce','Pr',
	                       'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
	                       'Yb','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk',
	                       'Cf','Es','Fm','Md','No'
	                       ]
	            for i in ele_symbols_list:
	                for j in metals:
	                    if j == i:
	                        return true

	            return false

	        add = 0
	        if contains_poly:
	            add = 1
	        for i in range(len(ele_symbols_list)+add):
	            ox_nums.append('0')
	        if contains_poly:
	            ox_nums.remove('0')
	            if poly_charge > 0:
	                ox_nums.append('+'+str(poly_charge))
	            else:
	                ox_nums.append(str(poly_charge))
	        for i in ele_symbols_list:
	            if i == 'O' and not 'F' in ele_symbols_list:
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '-2')
	            if i == 'F':
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '-1')
	            if (i == 'Cl' or i == 'Br' or i == 'I') and ('O' not in ele_symbols_list and 'F' not in ele_symbols_list):
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '-1')
	            if i == 'H' and not contains_metal(self):
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '+1')
	            if i == 'H' and contains_metal(self) and len(ele_symbols_list) == 2:
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '-1')
	            if i in alkaline_metals:
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '+1')
	            if i in alkaline_earth_metals:
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '+2')
	            if i == 'Al':
	                ox_nums.remove('0')
	                ox_nums.insert(ele_symbols_list.index(i), '+3')


	        vars_list = []
	        letters = 'abcdefghijklmnopqrstuvwxyz'
	        k = 0
	        index = 0
	        if len(digits) + 1 == len(ox_nums):
	            digits.append(None)
	        for i in ox_nums:
	            index = ox_nums.index(i)
	            factor = 0
	            if digits[index] == None:
	                factor = 1
	            else:
	                factor = digits[index]
	            if i == '0':
	                vars_list.append(str(factor)+'*'+letters[k])
	                k += 1
	        equation = ''
	        for i in vars_list:
	            equation += i + "+"
	        for i in ox_nums:
	            factor = 0
	            if digits[ox_nums.index(i)] == None:
	                factor = 1
	            else:
	                factor = digits[ox_nums.index(i)]
	            if i != '0':
	               equation += str(factor*int(i)) + "+"
	        equation = equation[:-1]
	        if has_charge and '0' in ox_nums:
	            if sign == '+':
	                equation += '-'+str(charge)
	            elif sign == '-':
	                equation += '+'+str(charge)
	            ox_nums[ox_nums.index('0')] = str(int(sympy.solve(equation, letters[k-1])[0]))
	        if '0' in ox_nums:
	            ox_nums[ox_nums.index('0')] = str(int(sympy.solve(equation, letters[k-1])[0]))
	        for i in range(len(ele_symbols_list)):
	            if i != len(ele_symbols_list) - 1:
	                self.result.text += ele_symbols_list[i] + ': ' + ox_nums[i] + ', '
	            else:
	                self.result.text += ele_symbols_list[i] + ': ' + ox_nums[i]

	        self.result.text = self.result.text.replace(', ', '')


	        for i in range(len(self.result.text)-1):
	            if self.result.text[i].isdigit() and self.result.text[i+1].isalpha():
	                replace = self.result.text[i] + self.result.text[i+1]
	                replacement = self.result.text[i] + ', ' +  self.result.text[i+1]
	                self.result.text = self.result.text.replace(replace, replacement)
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""

class PeriodicTable(Screen):
	pass
class LimitExcessReagent(Screen):
	pass
class Stoichiometry(Screen):
	equation = ObjectProperty(None)
	start = ObjectProperty(None)
	end = ObjectProperty(None)
	amount = ObjectProperty(None)
	start_unit = ObjectProperty(None)
	end_unit = ObjectProperty(None)
	result = ObjectProperty(None)

	item_index_a = 0
	item_index_b = 0

	units = ['grams', 'moles', 'liters', 'milliliters', 'particles']

	def change_item_a(self):
	    index = self.item_index_a
	    self.item_index_a += 1
	    if index > len(self.units) - 1:
	        index = 0
	        self.item_index_a = 1
	    self.start_unit.text = self.units[index]

	def change_item_b(self):
	    index = self.item_index_b
	    self.item_index_b += 1
	    if index > len(self.units) - 1:
	        index = 0
	        self.item_index_b = 1
	    self.end_unit.text = self.units[index]

	def read(self):
		try:
			start_comp = decode(self.start.text)
			end_comp = decode(self.end.text)
			start_u = self.start_unit.text
			end_u = self.end_unit.text
			amt = float(self.amount.text)
			self.start.text = ""
			self.end.text = ""
			self.amount.text = ""
			self.start_unit.text = 'Start Units'
			self.end_unit.text = 'End Units'
			eq = bal(decode(self.equation.text))[1]
			self.equation.text = ""
			def mass(self):
				start_mm = round(compounds.mass(start_comp), 1)
				ratio = eq[end_comp] / eq[start_comp]
				if end_u == 'moles':
					return str(amt / start_mm * ratio) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'liters' or end_u == "milliliters":
					if end_u == 'liters':
						return str(amt / start_mm * ratio * 22.4) + " " + end_u + " " + uncode(end_comp, 'subscripts')
					else:
						return str(amt / start_mm * ratio * 22400) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'grams':
					return str(amt / start_mm * ratio * round(compounds.mass(end_comp), 1)) + " " + end_u + "/moles" + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'particles':
					return str(amt / start_mm * ratio * (6.02 * 10**23)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
			def volume(self):
				ratio = eq[end_comp] / eq[start_comp]
				if end_u == 'moles':
					if start_u == 'liters':
						return str(amt / 22.4 * ratio) + " " + end_u + " " + uncode(end_comp, 'subscripts')
					else:
						return str(amt / 22400 * ratio) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'liters' or end_u == 'milliliters':
					if start_u == 'liters':
						if end_u == 'liters':
							return str(amt / 22.4 * ratio * 22.4) + " " + end_u + " " + uncode(end_comp, 'subscripts')
						else:
							return str(amt / 22.4 * ratio * 22400) + " " + end_u + " " + uncode(end_comp, 'subscripts')
					else:
						if end_u == 'liters':
							return str(amt / 22400 * ratio * 22.4) + " " + end_u + " " + uncode(end_comp, 'subscripts')
						else:
							return str(amt / 22400 * ratio * 22400) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'grams':
					if start_u == 'liters':
						return str(amt / 22.4 * ratio * round(compounds.mass(end_comp), 1)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
					else:
						return str(amt / 22400 * ratio * round(compounds.mass(end_comp), 1)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'particles':
					if start_u == 'liters':
						return str(amt / 22.4 * ratio * (6.02 * 10**23)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
					else:
						return str(amt / 22400 * ratio * (6.02 * 10 ** 23)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
			def particles(self):
				ratio = eq[end_comp] / eq[start_comp]
				if end_u == 'particles':
					return str(amt * ratio) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'grams':
					return str(amt / (6.02 * 10**23) * ratio * round(compounds.mass(end_comp), 1)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'liters' or end_u == 'milliliters':
					if end_u == 'liters':
						return str(amt / (6.02 * 10**23) * ratio * 22.4) + " " + end_u + " " + uncode(end_comp, 'subscripts')
					else:
						return str(amt / (6.02 * 10 ** 23) * ratio * 22400) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'moles':
					return str(amt / (6.02 * 10**23) * ratio) + " " + end_u + " " + uncode(end_comp, 'subscripts')
			def moles(self):
				ratio = eq[end_comp] / eq[start_comp]
				if end_u == 'moles':
					return str(amt * ratio) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'grams':
					return str(amt * ratio * round(compounds.mass(end_comp), 1)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'liters' or end_u == 'milliliters':
					if end_u == 'liters':
						return str(amt * ratio * 22.4) + " " + end_u + " " + uncode(end_comp, 'subscripts')
					else:
						return str(amt * ratio * 22400) + " " + end_u + " " + uncode(end_comp, 'subscripts')
				elif end_u == 'particles':
					return str(amt * ratio * (6.02 * 10**23)) + " " + end_u + " " + uncode(end_comp, 'subscripts')
			if start_u == 'grams':
				self.result.text = mass(self)
			elif start_u == 'liters' or start_u == 'milliliters':
				self.result.text = volume(self)
			elif start_u == 'particles':
				self.result.text = particles(self)
			elif start_u == 'moles':
				self.result.text = moles(self)
		except:
			self.result.text = "Invalid Input"
			self.equation.text = ""
			self.start.text = ""
			self.end.text = ""
			self.amount.text = ""
			self.start_unit.text = "Start Units"
			self.end_unit.text = "End Units"
	    
class Naming(Screen):
	type = ObjectProperty(None)
	data = ObjectProperty(None)
	nam = ObjectProperty(None)
	result = ObjectProperty(None)

	def assign(self):
	    try:
	        data_text = decode(self.data.text)
	        name_text = self.nam.text
	        type_text = self.type.text
	        self.data.text = ""
	        self.nam.text = ""
	        self.type.text = "molecular"
	        def ide(self, elements):
	            if elements[-1] == 'oxygen':
	                elements[-1] = 'oxide'
	            elif elements[-1] == 'hydrogen':
	                elements[-1] = 'hydride'
	            elif elements[-1] == 'nitrogen':
	                elements[-1] == 'nitride'
	            else:
	                s = elements[-1]
	                s = s[0:-3]
	                s += 'ide'
	                elements[-1] = s
	            return elements
	        def un_ide(self, elements):
	            s = elements[-1]
	            s = s[0:-3]
	            for i in ele_names:
	                if s in i:
	                    elements[-1] = i
	            return elements
	        def molecular(self):
	            prefixes = {'1':'mono','2':'di','3':'tri','4':'tetra','5':'penta',
	                        '6':'hexa','7':'hepta','8':'octa','9':'nona','10':'deca'
	                        }
	            completed = ''
	            if name_text.upper() == 'NAME':
	                pile = ""
	                for letter in data_text:
	                    pile = pile + letter + " "
	                pile = pile[:-1]
	                digits = list(int(s) for s in pile.split() if s.isdigit())
	                for i in digits:
	                    if i == 0:
	                        digits.remove(0)
	                        break
	                for i in range(len(digits)):
	                    if digits[i] == 1:
	                         digits[i] = 10
	                ele_symbols_str = ''.join(filter(str.isalpha, data_text))
	                ele_symbols_list = []
	                for i in range(len(ele_symbols_str)):
	                    if ele_symbols_str[i].isupper() and i != len(ele_symbols_str)-1 and ele_symbols_str[i+1].islower():
	                        ele_symbols_list.append(ele_symbols_str[i] + ele_symbols_str[i+1])
	                    elif ele_symbols_str[i].isupper():
	                        ele_symbols_list.append(ele_symbols_str[i])
	                if data_text.index(ele_symbols_list[0])+1 == data_text.index(ele_symbols_list[-1]):
	                    digits.insert(0, 1)
	                if data_text.index(ele_symbols_list[-1]) == len(data_text)-1:
	                    digits.append(1)
	                element_names = []
	                for i in range(len(ele_symbols_list)):
	                    for j in list(symbols.values()):
	                        if ele_symbols_list[i] == j:
	                            element_names.append(ele_names[symbols_list.index(j)])
	                element_names_revised = ide(self, element_names)
	                vowels = ['a', 'e', 'o', 'u']
	                list_pre = [prefixes[str(digits[0])], prefixes[str(digits[1])]]
	                for i in vowels:
	                    if i == element_names_revised[0][0]:
	                        list_pre[0] = list_pre[0][0:-1]
	                    if i == element_names_revised[1][0]:
	                        list_pre[1] = list_pre[1][0:-1]
	                if digits[0] == 1:
	                    completed += element_names_revised[0] + ' ' + list_pre[1] + element_names_revised[1]
	                else:
	                    completed += list_pre[0] + element_names_revised[0] + ' ' + list_pre[1] + element_names_revised[1]
	            elif name_text.upper() == 'FORMULA':
	                digit_list = [1, 1]
	                data_split = data_text.split()
	                prefix_keys = list(prefixes.keys())
	                prefix_values = list(prefixes.values())
	                for i in prefix_values:
	                    if i in data_split[0]:
	                        digit_list[0] = int(prefix_keys[prefix_values.index(i)])
	                        data_split[0] = data_split[0].replace(i, '')
	                    if i in data_split[1]:
	                        digit_list[1] = int(prefix_keys[prefix_values.index(i)])
	                        data_split[1] = data_split[1].replace(i, '')
	                symb_list = ['','']
	                element_list = un_ide(self, data_split)
	                symb_list[0] = symbols_list[ele_names.index(element_list[0])]
	                symb_list[1] = symbols_list[ele_names.index(element_list[1])]

	                if digit_list[0] == 1 and digit_list[1] == 1:
	                    completed = symb_list[0] + symb_list[1]
	                elif digit_list[0] == 1:
	                    completed = symb_list[0] + symb_list[1] + subscripts[str(digit_list[1])]
	                elif digit_list[1] == 1:
	                    completed = symb_list[0] + subscripts[str(digit_list[0])] + symb_list[1]
	                else:
	                    completed = symb_list[0] + subscripts[str(digit_list[0])] + symb_list[1] + subscripts[str(digit_list[1])]

	            return completed
	        def needsRomanNumerals(self, info, charges):
	            roman_numeral_elements = ['Sb','Bi','Au','Sn','Hg','Ni','Pb','Cu','Co','Fe','Mn','Cr'] # need more
	            for i in roman_numeral_elements:
	                if i == info[0]:
	                    info.append("("+int_to_roman(self,charges[info[0]])+")")
	                if i == info[1]:
	                    info.append("(" + int_to_roman(self, charges[info[1]]) + ")")
	            if len(info) > 2:
	                return true
	            else:
	                return false
	        def int_to_roman(self, input):
	            if not isinstance(input, type(1)):
	                raise TypeError("expected integer, got %s" % type(input))
	            if not 0 < input < 4000:
	                raise ValueError("Argument must be between 1 and 3999")
	            ints = (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1)
	            nums = ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
	            result = []
	            for i in range(len(ints)):
	                count = int(input / ints[i])
	                result.append(nums[i] * count)
	                input -= ints[i] * count
	            return ''.join(result)

	        def roman_to_int(self, input):
	            if not isinstance(input, type("")):
	                raise TypeError("expected string, got %s" % type(input))
	            input = input.upper(  )
	            nums = {'M':1000, 'D':500, 'C':100, 'L':50, 'X':10, 'V':5, 'I':1}
	            sum = 0
	            for i in range(len(input)):
	                try:
	                    value = nums[input[i]]
	                    # If the next place holds a larger number, this value is negative
	                    if i+1 < len(input) and nums[input[i+1]] > value:
	                        sum -= value
	                    else: sum += value
	                except KeyError:
	                    raise ValueError('input is not a valid Roman numeral: %s' % input)
	            # easiest test for validity...
	            if int_to_roman(sum) == input:
	                return sum
	            else:
	                raise ValueError('input is not a valid Roman numeral: %s' % input)
	        def getCharges(self, info):
	            temp = data_text
	            temp.replace(info[0],'')
	            temp.replace(info[1],'')
	            chrgs = list(int(s) for s in temp.split() if s.isdigit())
	            return chrgs
	        def ionic(self): # need to finish, lots of work to be done
	            completed = ''
	            if name_text.upper() == 'NAME':
	                info_list = []
	                for i in symbols_list:
	                    if i in data_text:
	                        info_list.append(i)
	                for j in polyatomics_keys:
	                    if j in data_text:
	                        info_list.append(j)
	                charges = getCharges(self, info_list)
	                if needsRomanNumerals(self, info_list, charges):
	                    completed += ele_names[symbols_list.index(info_list[0])] + ' ' + info_list[2] + ' ' + ele_names[symbols_list.index(info_list[1])]
	                else:
	                    completed += ele_names[symbols_list.index(info_list[0])] + ' ' + ele_names[symbols_list.index(info_list[1])]
	            elif name_text.upper() == 'FORMULA':
	                pass

	            return completed
	        if type_text.upper() == 'MOLECULAR':
	            self.result.text = molecular(self)
	        elif type_text.upper() == 'IONIC':
	            self.result.text = ionic(self)
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""
	        self.nam.text = ""

class Math(Screen):
	pass
class Calculus(Screen):
	pass
class Limits(Screen):
	expression = ObjectProperty(None)
	variable = ObjectProperty(None)
	approaches = ObjectProperty(None)
	result = ObjectProperty(None)

	def findTheLimit(self):
	    try:
	        ex = decode(self.expression.text, '^')
	        var = self.variable.text
	        apr = self.approaches.text
	        self.expression.text = ""
	        self.variable.text = ""
	        self.approaches.text = ""
	        self.result.text = str(limit(ex,var,apr))
	    except:
	        self.result.text = "Invalid Input"
	        self.expression.text = ""
	        self.variable.text = ""
	        self.approaches.text = ""
class Derivatives(Screen):
	pass
class SolveEquations(Screen):
	equation = ObjectProperty(None)
	variable = ObjectProperty(None)
	result = ObjectProperty(None)

	def solveeq(self):
	    try:
	        pi_symbol = u'\u03C0'
	        eq = self.equation.text
	        var = self.variable.text
	        eq = eq.replace(' ', '')
	        if '=0' in eq:
	            eq = eq.replace('=0', '')
	        elif '=0' not in eq and '=' in eq:
	            eq_list = eq.split('=')
	            newstr = ''
	            if eq_list[1][0] != '-':
	                eq_list[1] = eq_list[1].replace(eq_list[1][0], '+'+eq_list[1][0])
	            if '-' in eq_list[1]:
	                eq_list[1] = eq_list[1].replace('-', '+')
	            if '+' in eq_list[1]:
	                eq_list[1] = eq_list[1].replace('+','-')

	            eq = eq_list[0] + eq_list[1]

	        self.equation.text = ""
	        self.variable.text = ""
	        eq = decode(eq, '^')
	        for i in range(len(eq)-1):
	            if eq[i].isdigit() and eq[i+1].isalpha():
	                replacement = eq[i] + "*" + eq[i+1]
	                replaced = eq[i] + eq[i+1]
	                eq = eq.replace(replaced, replacement)
	        solved = str(sympy.solve(eq, var)).replace('sqrt', sqrt_symbol)
	        solved = solved.replace('pi', pi_symbol)
	        self.result.text = solved
	    except:
	        self.result.text = "Invalid Input"
	        self.equation.text = ""
	        self.variable.text = ""
class Trigonometry(Screen):
	pass
	degrees = pi/180
class UnitCircle(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def evaulateTrig(self):
	    try:
	        data_list = self.data.text.split(',')
	        self.data.text = ""

	        func = ''
	        ang = 0
	        if 'sin' in data_list[0]:
	            func = 'sin'
	        elif 'cos' in data_list[0]:
	           func = 'cos'
	        elif 'tan' in data_list[0]:
	            func = 'tan'
	        elif 'cot' in data_list[0]:
	            func = 'cot'
	        elif 'csc' in data_list[0]:
	            func = 'csc'
	        elif 'sec' in data_list[0]:
	            func = 'sec'

	        newstr = ''.join((ch if ch in 'piπ0123456789/.-' else ' ') for ch in data_list[1])
	        if 'pi' in newstr:
	            newstr = newstr.replace('pi', ' pi ')
	        elif 'π' in newstr:
	            newstr = newstr.replace('π', ' π ')
	        elif '/' in newstr:
	            newstr = newstr.replace('/', ' / ')
	        temp = newstr.split()
	        if len(temp) > 1:
	            if temp[0].isdigit() and temp[1] == 'π' or temp[1] == 'pi':
	                data_list[1] = data_list[1].replace(temp[0]+temp[1], temp[0]+'*'+temp[1])


	        if 'angle' in data_list[1] and '°' in data_list[1]:
	            data_list[1] = data_list[1].replace('°', '*'+str(degrees))
	        elif 'angle' in data_list[1]:
	            if 'π' in data_list[1]:
	                data_list[1] = data_list[1].replace('π', str(pi))
	            elif 'pi' in data_list[1]:
	                data_list[1] = data_list[1].replace('pi', str(pi))

	        i = data_list[1]
	        if '*' in i and '/' in i and i.index('*') < i.index('/'):
	            ang = float(i[i.index('=')+1:i.index('*')]) * float(i[i.index('*')+1:i.index('/')]) / float(i[i.index('/')+1:])
	        elif '*' in i and '/' in i and i.index('*') > i.index('/'):
	            ang = float(i[i.index('=')+1:i.index('/')]) / float(i[i.index('/')+1:i.index('*')]) * float(i[i.index('*')+1:])
	        elif '/' in i and not '*' in i:
	            ang = float(i[i.index('=')+1:i.index('/')]) / float(i[i.index('/')+1:])
	        elif '*' in i and not '/' in i:
	            ang = float(i[i.index('=')+1:i.index('*')]) * float(i[i.index('*')+1:])

	        if func == 'sin':
	            self.result.text = str(sympy.simplify(sympy.sin(ang)))
	        elif func == 'cos':
	            self.result.text = str(sympy.simplify(sympy.cos(ang)))
	        elif func == 'tan':
	            self.result.text = str(sympy.simplify(sympy.tan(ang)))
	        elif func == 'csc':
	            self.result.text = str(sympy.simplfy(1/sympy.sin(ang)))
	        elif func == 'sec':
	            self.result.text = str(sympy.simplify(1/sympy.cos(ang)))
	        elif func == 'cot':
	            self.result.text = str(sympy.simplify(1/sympy.tan(ang)))
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""

class RightTriangleTrig(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def solve(self):
	    try:
	        data_text = self.data.text
	        self.data.text = ""

	        newstr = ''.join((ch if ch in '°*degreespiπ/0123456789.-e' else ' ') for ch in data_text)

	        for i in newstr.split():
	            if 'pi' in i:
	                newstr = newstr.replace('pi',str(pi))
	            if 'π' in i:
	                newstr = newstr.replace('π', str(pi))
	            if 'degrees' in i:
	                newstr = newstr.replace('degrees', str(degrees))
	            if '°' in i:
	                newstr = newstr.replace('°', '*'+str(degrees))

	        div_list = []
	        mult_list = []
	        nums_div = []
	        nums_mult = []
	        temp = newstr.split()
	        for i in temp:
	            if '/' in i:
	                div_list.append(i)
	                if str(pi) in i:
	                    nums_div += temp[temp.index(i)].split('/')
	            if '*' in i:
	                mult_list.append(i)
	                if str(degrees) in i:
	                    nums_mult += temp[temp.index(i)].split('*')


	        j = 0
	        for i in range(len(div_list)):
	            newstr = newstr.replace(div_list[i], str(float(nums_div[j]) / float(nums_div[j+1])))
	            j += 2

	        k = 0
	        for i in range(len(mult_list)):
	            newstr = newstr.replace(mult_list[i], str(float(nums_mult[k]) * float(nums_mult[k+1])))
	            k += 2

	        vals = [float(i) for i in newstr.split()]

	        letters = []
	        for character in data_text:
	            if character.isalpha() and character not in 'pi' and character not in 'degrees' and character != 'π':
	                letters.append(character)

	        lets_list = ['a', 'b', 'c', 'A', 'B', 'C']
	        args_list = [None, None, None, None, None, None]

	        for i in letters:
	            for j in lets_list:
	                if i == j:
	                    args_list.insert(lets_list.index(j), vals[letters.index(i)])

	        ans_tup = solve(a=args_list[0], b=args_list[1], c=args_list[2], A=args_list[3], B=args_list[4], C=args_list[5])

	        self.result.text = 'a = '+str(ans_tup[0])+'\nb = '+str(ans_tup[1])+'\nc = '+str(ans_tup[2])+'\nA = '+str(ans_tup[3])+'\nB = '+str(ans_tup[4])+'\nC = '+str(ans_tup[5])
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""

class Factor(Screen):
	data = ObjectProperty(None)
	result = ObjectProperty(None)

	def factorize(self):
	    try:
	        data_text = decode(self.data.text, '^')
	        self.data.text = ""
	        for i in range(len(data_text)-1):
	            exp = data_text[i]+data_text[i+1]
	            if data_text[i].isdigit() and data_text[i+1].isalpha():
	                newstr = data_text[i]+'*'+data_text[i+1]
	                data_text = data_text.replace(exp, newstr)
	        ans = str(sympy.factor(data_text))
	        for i in range(len(ans)-2):
	            if ans[i] == '*' and ans[i+1] == '*':
	                digit = ans[i+2]
	                ans = ans.replace(ans[i:i+3], superscripts[digit])
	        self.result.text = ans.replace('*','')
	    except:
	        self.result.text = "Invalid Input"
	        self.data.text = ""
class Balance(Screen):
	equation = ObjectProperty(None)
	result = ObjectProperty(None)

	def balance(self):
		try:
			self.result.text = bal(decode(self.equation.text))[0]
			self.equation.text = ""
		except:
			self.result.text = "Invalid Input"
			self.equation.text = ""

class RadicalOperations(Screen):
	pass
class MainApp(MDApp):
	dialog1 = None
	dialog2 = None
	dialog3 = None
	dialog4 = None
	dialog5 = None
	dialog6 = None
	dialog7 = None
	dialog8 = None
	dialog9 = None
	dialog10 = None
	dialog11 = None
	dialog12 = None
	dialog13 = None
	dialog14 = None
	dialog15 = None
	dialog16 = None

	def build(self):
	    self.theme_cls.primary_palette = "BlueGray"
	    self.theme_cls.theme_style = "Dark"
	    self.theme_cls.accent_hue = "900"

	    return Builder.load_file('main.kv')

	def show_theme_picker(self):
	    theme_dialog = MDThemePicker()
	    theme_dialog.open()

	def show_MDDialog_1(self, text):
	    if not self.dialog1:
	        self.dialog1 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog1.open()

	def show_MDDialog_2(self, text):
	    if not self.dialog2:
	        self.dialog2 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog2.open()

	def show_MDDialog_3(self, text):
	    if not self.dialog3:
	        self.dialog3 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog3.open()

	def show_MDDialog_4(self, text):
	    if not self.dialog4:
	        self.dialog4 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog4.open()

	def show_MDDialog_5(self, text):
	    if not self.dialog5:
	        self.dialog5 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog5.open()

	def show_MDDialog_6(self, text):
	    if not self.dialog6:
	        self.dialog6 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog6.open()

	def show_MDDialog_7(self, text):
	    if not self.dialog7:
	        self.dialog7 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog7.open()

	def show_MDDialog_8(self, text):
	    if not self.dialog8:
	        self.dialog8 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog8.open()

	def show_MDDialog_9(self, text):
	    if not self.dialog9:
	        self.dialog9 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog9.open()

	def show_MDDialog_10(self, text):
	    if not self.dialog10:
	        self.dialog10 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog10.open()

	def show_MDDialog_11(self, text):
	    if not self.dialog11:
	        self.dialog11 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog11.open()

	def show_MDDialog_12(self, text):
	    if not self.dialog12:
	        self.dialog12 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog12.open()

	def show_MDDialog_13(self, text):
	    if not self.dialog13:
	        self.dialog13 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog13.open()

	def show_MDDialog_14(self, text):
	    if not self.dialog14:
	        self.dialog14 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog14.open()

	def show_MDDialog_15(self, text):
	    if not self.dialog15:
	        self.dialog15 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog15.open()

	def show_MDDialog_16(self, text):
	    if not self.dialog16:
	        self.dialog16 = MDDialog(
	            title="Info",
	            text=text,
	            size_hint=[0.5, 0.5],
	            pos_hint={'center_x': 0.5, 'center_y': 0.5},
	        )
	    self.dialog16.open()


if __name__ == "__main__":
	MainApp().run()