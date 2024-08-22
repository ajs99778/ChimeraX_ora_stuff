from chimerax.atomic import AtomicStructure
from chimerax.core.commands import run
import numpy as np

def run_preset(session, name, mgr):
    if name == "ORA style":
        ora_style(session)
    else:
        raise ValueError

def ora_style(session):
    # iORA's default lighting color is darkgray
    run(session, "lighting color darkgray ambientColor darkgray fillColor darkgray depthCue false", log=False)
    # iORA's atom colors are pretty close to Jmol's (aka byatom coloring in ChimeraX)
    run(session, "color byatom", log=False)
    # except carbons are black
    run(session, "color C black", log=False)
    run(session, "material reflectivity 0.75", log=False)
    
    radii = {
        "Ag": 0.651215197977908,
        "Al": 0.5529743063850676,
        "Ar": 0.48516474400095355,
        "As": 0.5560183626608767,
        "Au": 0.6274603577123187,
        "B": 0.4315147710181497,
        "Ba": 0.7575567729333123,
        "Be": 0.46072857260889966,
        "Bi": 0.6328217606804789,
        "Br": 0.5406399435280995,
        "C": 0.41250000000000003,
        "Ca": 0.7032122679680467,
        "Cd": 0.638135293340381,
        "Cl": 0.4919660946512383,
        "Co": 0.5769016178022384,
        "Cr": 0.5798260055265951,
        "Cs": 0.8132925069349483,
        "Cu": 0.6110801024190056,
        "F": 0.388850065285295,
        "Fe": 0.5739627885984686,
        "Ga": 0.5769016178022384,
        "Ge": 0.5650582079376084,
        "H": 0.23351768280630916,
        "He": 0.20682658593383532,
        "Hf": 0.6434018029708556,
        "Hg": 0.6407743739596182,
        "I": 0.5970768280006289,
        "In": 0.6274603577123187,
        "Ir": 0.6083056838538263,
        "K": 0.7532128725927925,
        "Kr": 0.5280445166670034,
        "La": 0.6912377854157897,
        "Li": 0.5999039673050671,
        "Lu": 0.6690571421338852,
        "Mg": 0.5885139219318248,
        "Mn": 0.613841646605784,
        "Mo": 0.6301470968327719,
        "N": 0.4047206512588323,
        "Na": 0.6537970315474769,
        "Nb": 0.6083056838538263,
        "Ne": 0.3807531825013603,
        "Ni": 0.5620601577461922,
        "O": 0.3968382617929041,
        "Os": 0.5827360930093497,
        "P": 0.5151767352555208,
        "Pb": 0.6354844573121149,
        "Pd": 0.5913819356965537,
        "Pt": 0.5827360930093497,
        "Rb": 0.785044615058307,
        "Re": 0.6665408098251128,
        "Rh": 0.602717739470673,
        "Rn": 0.6301470968327719,
        "Ru": 0.5769016178022384,
        "S": 0.5020245600848591,
        "Sb": 0.6110801024190056,
        "Sc": 0.6274603577123187,
        "Se": 0.5468390798364364,
        "Si": 0.5312184532803786,
        "Sn": 0.6193265859338354,
        "Sr": 0.7444287765179607,
        "Ta": 0.6110801024190056,
        "Tc": 0.6589273366332439,
        "Te": 0.602717739470673,
        "Ti": 0.6055182703065116,
        "Tl": 0.638135293340381,
        "V": 0.5739627885984686,
        "W": 0.6328217606804789,
        "Xe": 0.5885139219318248,
        "Y": 0.6740581102554121,
        "Zn": 0.5913819356965537,
        "Zr": 0.638135293340381
    }
    
    
    for m in session.models.list(type=AtomicStructure):
        try:
            m.ball_scale = 1
            for bond in m.bonds:
                bond.halfbond = False
                bond.color = np.array([255, 255, 255, 255], dtype=np.uint8)
                bond.radius = 0.06
            
            for atom in m.atoms:
                atom.draw_mode = atom.BALL_STYLE
                try:
                    atom.radius = radii[atom.element.name]
                except KeyError:
                    session.logger.warning("no radius for", atom.element.name)
                    atom.radius = radii["H"]

        except AttributeError:
            print("cannot color bond")
    