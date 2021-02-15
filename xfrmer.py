import math as m
from collections import namedtuple

EI = namedtuple('EI', ['name', 'core_width', 'winding_area_width', 'winding_area_height', 'max_height'])

EI_cores = [
    EI("EI35", 0.0096, 0.0195, 0.0077,0.0192),
    EI("EI42", 0.014, 0.021, 0.008,0.028),
    EI("EI54", 0.018, 0.027, 0.009,0.036),
    EI("EI60", 0.020, 0.030, 0.01,0.04),
    EI("EI76.2A", 0.0254, 0.0381, 0.0127,0.0508),
    EI("EI96", 0.032, 0.064, 0.048,0.064),
    EI("EI105", 0.035, 0.0525, 0.0175,0.07),
    EI("EI114.3", 0.0381, 0.05715, 0.01905,0.0762),
    EI("EI133.2", 0.04445, 0.0666, 0.0222,0.0889)    
    ]

Wire = namedtuple('Wire', ['name', 'cond_diameter', 'total_diameter'])
Wires = [
    Wire("10#",0.002588,0.002703),
    Wire("11#",0.002304,0.002418),
    Wire("12#",0.002052,0.002163),
    Wire("13#",0.001829,0.001934),
    Wire("14#",0.001628,0.001732),
    Wire("15#",0.00145,0.001547),
    Wire("16#",0.00129,0.001384),
    Wire("17#",0.001151,0.00124),
    Wire("18#",0.001024,0.00111),
    Wire("19#",0.000912,0.000993),
    Wire("20#",0.000813,0.000892),
    Wire("21#",0.000724,0.000798),
    Wire("22#",0.000643,0.000714),
    Wire("23#",0.000574,0.000642),
    Wire("24#",0.000511,0.000577),
    Wire("25#",0.000455,0.000516),
    Wire("26#",0.000404,0.000462),
    Wire("27#",0.000361,0.000417),
    Wire("28#",0.00032,0.000373),
    Wire("29#",0.000287,0.000338),
    Wire("30#",0.000254,0.000302),
    Wire("31#",0.000226,0.000274),
    Wire("32#",0.000203,0.000249),
    Wire("33#",0.00018,0.000224),
    Wire("34#",0.00016,0.000198),
    Wire("35#",0.000142,0.000178),
    Wire("36#",0.000127,0.00016),
    Wire("37#",0.000114,0.000145),
    Wire("38#",0.000102,0.00013),
    Wire("39#",0.000089,0.000114),
    Wire("40#",0.000079,0.000102),
    Wire("41#",0.000071,0.000091),
    Wire("42#",0.000064,0.000081),
    Wire("43#",0.000056,0.000074),
    Wire("44#",0.000051,0.000069),
    Wire("45#",0.0000447,0.000061),
    Wire("46#",0.0000399,0.0000533)
    ]
Soln = namedtuple('Soln', ['wp', 'ws', 'core', 'stack_height'])
Solns = []
CU_rho = 1.68 * 10**(-8)

def Afromd(diameter):
    return m.pi*(diameter/2)**2

def turns_per_layer(wire_size, width):
    return m.floor(width/wire_size)

def number_of_layers(tpl, turns):
    return m.ceil(turns/tpl) 
   
def length_of_wire(turns, stack_height, core_width, winding_area_width, wire_size):
    rt = 0 #running total
    tpl = turns_per_layer(wire_size, winding_area_width)
    while(turns > 0):
        if(turns > tpl):
            rt += tpl*2*(core_width + stack_height)
            turns -= tpl
            core_width += 2*wire_size
            stack_height += 2*wire_size
        elif(turns <= tpl and turns > 0):
            rt += 2*turns*(core_width + stack_height)
            turns = 0
    return rt

def inputs():
    J_in = 0 #input Current density [A/mm^2]
    J_out = 0 #output Current density [A/mm^2]
    B_sat = 0 #Magnetic saturation point [T]
    I_in_max = 0 #Maximum input current [A]
    I_out_max = 0 #Maximum output current [A]
    V_in = 0 #Input voltage
    outputs = 0 #number of outputs
    V_out = [] #Output voltages
    f = 0 #input frequency

def offset_height(r):
    return r*(1+m.sqrt(3))

def nturns(V, area, f, Bsat):
    return V/4.44/area/f/Bsat

def zipmap(f, x):
    return zip(x, map(f, x))

def repacktuple(t, n):
    return list(zip(*zip(*t), n))

#def letsdosomeshit(I_out_max, ):
I_out_max = 5
V_out = 30
P_out = I_out_max * V_out
P_loss_max = P_out / 10

core = EI_cores[1]
solved = False

foundSoln = False
#Here we have made the assumption n2/n1 = v2/v1
R_max_p = P_loss_max/2/(V_out/220*I_out_max)**2
R_max_s = P_loss_max/2/I_out_max**2
#search space is quite small, we are iterating through available cores, looking for the smallest suitable wires to produce results.
for core in EI_cores:
    stack_height = 0.035
    foundSoln = False
    while(stack_height <= core.max_height and foundSoln == False):
        A = stack_height*core.core_width
        n2=nturns(V_out, A, 50, 1.2) #secondary side
        n1=nturns(220, A, 50, 1.2) #primary

        #calculate minimum wire diameter, this is assumeing output current > input current ie stepdown.
        
        #create a list of lengths for the differnt wire sizes
        lengths = list(map(lambda x: length_of_wire(n1, stack_height, core.core_width, core.winding_area_width, x.total_diameter), Wires))
        resistances = list(map(lambda x,y:  y/Afromd(x.cond_diameter)*CU_rho, Wires, lengths))
        smallest_wires = list(filter(lambda x: x[2] < R_max_p, zip(Wires, lengths, resistances)))
        smallest_wires.reverse()
        if(len(smallest_wires) >= 1):
            smallest_wire = smallest_wires[0][0]
        else:
            stack_height += 0.0005
            continue
        #height_of_wires = list(map(lambda x: x[0].total_diameter*number_of_layers(turns_per_layer(x[0].total_diameter, core.winding_area_height), n1), suitable_wires))
        height_of_wire =  smallest_wire.total_diameter*number_of_layers(turns_per_layer(smallest_wire.total_diameter, core.winding_area_height), n1)
        #First winding will fit in the core, check the second winding
        if(height_of_wire < core.winding_area_height):
            lengths2 = list(map(lambda x: length_of_wire(n2, stack_height+height_of_wire, core.core_width+height_of_wire, core.winding_area_width, x.total_diameter), Wires))
            resistances2 = list(map(lambda x,y:  y/Afromd(x.cond_diameter)*CU_rho, Wires, lengths2))
            smallest_wires2 = list(filter(lambda x: x[2] < R_max_s, zip(Wires, lengths2, resistances2)))
            smallest_wires2.reverse()
            if(len(smallest_wires2) >= 1):
                smallest_wire2 = smallest_wires2[0][0]
            else:
                stack_height += 0.0005
                continue
            height_of_wire2 =  smallest_wire2.total_diameter*number_of_layers(turns_per_layer(smallest_wire2.total_diameter, core.winding_area_height), n2)
            #Both windings will fit, assuming maximum fill of 90%
            if(height_of_wire + height_of_wire2 < core.winding_area_height*.9):
                Solns.append(Soln(smallest_wire, smallest_wire2, core, stack_height))
                foundSoln = True
        #add extra lamination
        stack_height += 0.0005
Solns


#suitable_wires = [w for w in resistances if w[2] < R_max]         

#check power loss
#length = length_of_wire(n1, stack_height, core.core_width, core.winding_area_width, wire.total_diameter)

#min_diameter = m.sqrt((4*length*CU_rho)/(R_max * m.pi))
#find smallest wire which 

    
        
    
#letsdosomeshit(5, 30)
#Questions:
    #should wires be wound straight on top of each other or offset?
    #Assuming widing primary first